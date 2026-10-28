// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>
// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using GlycoNetworks;
using MathNet.Numerics.Statistics;
using System.Diagnostics;

const double defaultMassTolerance = 0.03;
const int defaultVerbosity = 1;

if (args.Length < 2 || args.Length > 6)
{
    Console.WriteLine("Usage: GlycoNetworksDelta <Spectra file path> <Glycan database file path> [Mass tolerance] [Deltas]");
    Console.WriteLine("  <Spectra file path> is mandatory. Must be in MGF format, or ConvertToMgf.bat must be configured to convert to MGF format");
    Console.WriteLine("  <Glycan database file path> is mandatory. Must be a text file of glycan compositions");
    Console.WriteLine("  [Verbosity] is optional. Current verbosity levels are 0, 1, or 2. Level 0 outputs delta mass and confidence interval. Level 1 provides basic detail. Level 2 shows 'all' details. Default is: " + defaultVerbosity);
    Console.WriteLine("  [Mass tolerance] is optional. Default is: " + defaultMassTolerance);
    Console.WriteLine("  [Deltas] is optional. List of glycan deltas to consider, separated by commas. You may also specify a delta mass");
    Console.WriteLine("    Default is: all building blocks found in the glycan database, along with the disaccharide Hex-HexNAc");
    Console.WriteLine("    This program knows the following deltas: " + String.Join(",", Constants.glycan_masses.Keys));
    return;
}

string spectraFilePath = Path.GetFullPath(args[0]);

if (Path.GetExtension(spectraFilePath).ToLower() != ".mgf")
{
    var processStartInfo = new ProcessStartInfo(Path.Combine(AppContext.BaseDirectory, "ConvertToMgf.bat"), '"' + spectraFilePath + '"');
    processStartInfo.UseShellExecute = false;
    processStartInfo.CreateNoWindow = true;
    processStartInfo.WorkingDirectory = Path.GetDirectoryName(spectraFilePath);
    Process? process = Process.Start(processStartInfo);
    process.WaitForExit();

    spectraFilePath = Path.Combine(Path.GetDirectoryName(spectraFilePath), Path.GetFileNameWithoutExtension(spectraFilePath) + ".mgf");
    if (!File.Exists(spectraFilePath))
    {
        Console.WriteLine("Error: Unable to convert " + args[0] + " to MGF format");
        return;
    }
}

string glycanDatabaseFilePath = args[1];

if (args.Length >= 3)
{
    if (!int.TryParse(args[2], out Constants.verbosityLevel))
    {
        Console.WriteLine("Error: Could not parse weight");
        return;
    }
}
else
{
    Constants.verbosityLevel = defaultVerbosity;
}

if (args.Length >= 4)
{
    if (!Double.TryParse(args[3], out Constants.tolerance))
    {
        Console.WriteLine("Error: Could not parse mass tolerance");
        return;
    }
}
else
{
    Constants.tolerance = defaultMassTolerance;
}

string[]? deltas = null;
if (args.Length >= 5)
{
    deltas = args[4].Split(",");
    foreach (string delta in deltas)
    {
        if (!Constants.glycan_masses.ContainsKey(delta))
        {
            // In addition to pre-defined deltas, such as Hex or HexNAc, allow deltas of arbitrary mass (user must enter as floating point number)
            double num;
            if (Double.TryParse(delta, out num) && num != 0)
            {
                Constants.glycan_masses.Add(delta, num);
            }
            else
            {
                Console.WriteLine("Error: Invalid delta");
                return;
            }
        }
    }
}

string spectraFileName = Path.GetFileName(spectraFilePath);
string spectraDirectory = Path.GetDirectoryName(spectraFilePath);

Directory.CreateDirectory(Path.Combine(spectraDirectory, ("spectra - " + spectraFileName)));
Directory.CreateDirectory(Path.Combine(spectraDirectory, ("graphs - " + spectraFileName)));

List<GlycanComposition> compositions = Utilities.readGlycanDatabase(glycanDatabaseFilePath);
if (deltas == null)
{
    var deltasHashSet = new HashSet<string>();
    foreach (GlycanComposition composition in compositions)
    {
        foreach (string delta in composition.Keys)
        {
            deltasHashSet.Add(delta);
        }
    }
    if (deltasHashSet.Contains("Hex") && deltasHashSet.Contains("HexNAc"))
    {
        deltasHashSet.Add("Hex-HexNAc");
    }
    deltas = deltasHashSet.ToArray();
}
Console.WriteLine();
Console.WriteLine("Deltas: " + String.Join(",", deltas));
Console.WriteLine();

List<MS2> spectra = Utilities.readMgfAndFindGlycopeptides(Path.Combine(spectraDirectory, spectraFileName)); // reads in mgf, returns set of spectra (based on glycopeps)
var glycanMasses = new List<double>();
var glycopepList = new List<Glycopep>();  // Retention time, Charge, Precursor Mz, Precursor Mass, Peptide Mass, Glycan Mass, Intensity, [Hexes, HexNacs, Fucs, NeuAcs, NeuGcs]*
var commonGlycans = new List<GlycInfo>();

var mapGlycantoPepList = new Dictionary<double, List<Glycopep>>(); //glycmass, list of pepmasses
var mapPeptoGlycanList = new Dictionary<double, List<Glycopep>>();

var uniqueGlycansGlycopepList = new List<Glycopep>();

var otherGraph = new Dictionary<Glycopep, List<Edge<Glycopep>>>();
var graph3 = new Dictionary<Glycopep, List<Edge<Glycopep>>>();

using (var gnuplotConvenienceFile = new StreamWriter(Path.Combine(spectraDirectory, ("spectra - " + spectraFileName + @"\all.txt"))))
{
    using (var spectraSummaryCsvFile = new StreamWriter(Path.Combine(spectraDirectory, (spectraFileName + " - spectra summary.csv"))))
    {
        spectraSummaryCsvFile.WriteLine("Retention Time (min),Charge,Precursor Mz,Precursor Mass,Peptide Mass,Glycan Mass,Intensity");

        for (int i = 0; i < spectra.Count; i++)
        {
            Point? pepMz = spectra[i].findGlycopep();  // Look for pep, pep+HexNAc, pep+2HexNAc, etc. signature
            if (pepMz != null && pepMz.mass > 0)
            {
                double pepMass = Utilities.mzToMass(pepMz.mass, spectra[i].returnBarePepCharge());
                double glycMass = Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge) - pepMass;
                glycanMasses.Add(glycMass);

                spectraSummaryCsvFile.WriteLine((spectra[i].Rt / 60).ToString("F4") + "," + spectra[i].Charge + "," + spectra[i].PrecursorMz + ","
                    + Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge) + "," + pepMass + ","
                    + glycMass + "," + pepMz.intensity);

                Glycopep glycopep = new Glycopep();
                glycopep.rt = Convert.ToDouble((spectra[i].Rt / 60).ToString("F4"));
                glycopep.charge = spectra[i].Charge;
                glycopep.precursorMz = spectra[i].PrecursorMz;
                glycopep.precursorMass = Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge);
                glycopep.peptideMass = pepMass;
                glycopep.glycanMass = glycMass;
                glycopep.intensity = pepMz.intensity;
                glycopepList.Add(glycopep);

                string gnuplotSpectrumFileName = Path.Combine(spectraDirectory, @"spectra - " + spectraFileName + @"\rt=" + (spectra[i].Rt / 60).ToString("F4") + "_mz=" + spectra[i].PrecursorMz.ToString("F4") + ".txt");
                gnuplotConvenienceFile.WriteLine("load '" + Path.GetFileName(gnuplotSpectrumFileName) + "'; pause -1");

                using (var gnuplotSpectrumFile = new StreamWriter(gnuplotSpectrumFileName))
                {
                    string plotTitle = @"Retention time: " + spectra[i].Rt / 60 + @"\nCharge: " + spectra[i].Charge + @"\nPrecursor m/z: " + spectra[i].PrecursorMz;
                    gnuplotSpectrumFile.WriteLine("unset key\nset title \"" + plotTitle + "\"\nplot '-' with impulse linewidth 2 linetype - 1, '-' with impulse linewidth 2 linetype 7, '-' with impulse linewidth 2 linetype 17");

                    // Draw spectrum in black in gnuplot
                    for (int j = 0; j < spectra[i].findLen(); j++)
                    {
                        gnuplotSpectrumFile.WriteLine((spectra[i].GetMass(j)).ToString("F4") + " " + (spectra[i].GetIntensity(j)).ToString("F4"));
                    }

                    spectra[i].sortPeaksByIntensity();  

                    // Make the pep, pep+HexNAc, pep+2HexNAc, pep+2HexNAc+Hex peaks red in gnuplot
                    Point plusHexnac = spectra[i].returnHexnac();
                    Point plus2Hexnac = spectra[i].return2Hexnac();
                    Point plusHex2Hexnac = spectra[i].returnHex2Hexnac();
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n" + pepMz.mass.ToString("F4") + " " + pepMz.intensity);
                    gnuplotSpectrumFile.WriteLine(plusHexnac.mass + " " + plusHexnac.intensity);
                    gnuplotSpectrumFile.WriteLine(plus2Hexnac.mass + " " + plus2Hexnac.intensity);
                    gnuplotSpectrumFile.WriteLine(plusHex2Hexnac.mass + " " + plusHex2Hexnac.intensity);
                    gnuplotSpectrumFile.WriteLine("# " + (spectra[i].Rt / 60).ToString("F4"));

                    // Make the 204 and 366 peaks purple in gnuplot
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n" + Constants.glycan_masses["HexNAc"] + Constants.adduct + " " + spectra[i].corrIntensity(Constants.glycan_masses["HexNAc"] + Constants.adduct));
                    gnuplotSpectrumFile.WriteLine(Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct + " " + spectra[i].corrIntensity(Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct));
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n");
                }
            }
        }
    }
}

CompareGlycan glycans = new CompareGlycan();

for (int i = 0; i < glycopepList.Count; i++)
{
    for (int j = 0; j < compositions.Count; j++)
    {
        if (Math.Abs(glycopepList[i].glycanMass - glycans.calcMass(compositions[j])) <= Constants.tolerance)
        {
            glycopepList[i].comp = compositions[j];
        }
    }
}

for (int i = 0; i < glycopepList.Count; i++)
{
    if (Array.Exists(mapGlycantoPepList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].glycanMass) <= Constants.tolerance))
    {
        double keyVal = Array.Find(mapGlycantoPepList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].glycanMass) <= Constants.tolerance);
        List<Glycopep> list = mapGlycantoPepList[keyVal];
        list.Add(glycopepList[i]);
        mapGlycantoPepList[keyVal] = list;
    } else
    {
        List<Glycopep> list = new List<Glycopep>(); 
        list.Add(glycopepList[i]);
        mapGlycantoPepList.Add(glycopepList[i].glycanMass, list);
    }
 
    if (Array.Exists(mapPeptoGlycanList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].peptideMass) <= Constants.tolerance))
    {
        double keyVal = Array.Find(mapPeptoGlycanList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].peptideMass) <= Constants.tolerance);
        List<Glycopep> list = mapPeptoGlycanList[keyVal];
        list.Add(glycopepList[i]);
        mapPeptoGlycanList[keyVal] = list;       
    }
    else
    {
        List<Glycopep> list = new List<Glycopep>();
        list.Add(glycopepList[i]);
        mapPeptoGlycanList.Add(glycopepList[i].peptideMass, list);
    }

    GlycInfo temp = new GlycInfo();
    temp.glycMass = glycopepList[i].glycanMass;
    temp.occurrences = 1;
    temp.rtOcc = glycopepList[i].rt;
    temp.comp = glycopepList[i].comp;
    temp.clusterSize = 0;
    commonGlycans.Add(temp);

    for (int j = i + 1; j < glycopepList.Count; j++)
    {
        if ((Math.Abs(glycopepList[i].glycanMass - glycopepList[j].glycanMass) <= Constants.tolerance) && (Math.Abs(glycopepList[i].peptideMass - glycopepList[j].peptideMass) <= Constants.tolerance))
        {
            commonGlycans[i].occurrences++;
            if (glycopepList[j].comp != null)
            {
                commonGlycans[i].comp = glycopepList[j].comp;
            }
            glycopepList.RemoveAt(j);
            j--;
        }
    }
}

for (int i = 0; i < mapPeptoGlycanList.Count; i++)
{
    double keyVal = mapPeptoGlycanList.Keys.ToArray()[i];
    for (int k = 0; k < mapPeptoGlycanList[keyVal].Count; k++)
    {
        //mapPeptoGlycanList[keyVal].comps.Add("");
        for (int j = 0; j < glycopepList.Count; j++)
        {
            if (Math.Abs(mapPeptoGlycanList[keyVal][k].glycanMass - glycopepList[j].glycanMass) <= Constants.tolerance && glycopepList[j].comp != null)
            {
                mapPeptoGlycanList[keyVal][k].comp = (glycopepList[j].comp);
            }
        }
    }
}

for (int i = 0; i < mapGlycantoPepList.Count; i++)
{
    double keyVal = mapGlycantoPepList.Keys.ToArray()[i];
    uniqueGlycansGlycopepList.Add(mapGlycantoPepList[keyVal][0]);
}

// All-in-one graph is deprecated
using (var glycanGraphFile = new StreamWriter(Path.Combine(spectraDirectory, @"glycan graph - " + spectraFileName + ".gml")))
{
    glycanGraphFile.WriteLine("Creator \"yFiles\"\nVersion 2.2\ngraph\n[ hierarchic  1\n  directed  1  ");

    for (int i = 0; i < uniqueGlycansGlycopepList.Count; i++)
    {
        string label = "";
        if (uniqueGlycansGlycopepList[i].comp != null)
        {
            label = uniqueGlycansGlycopepList[i].comp.toString() + ": ";
        }
        glycanGraphFile.WriteLine("  node");
        glycanGraphFile.WriteLine("  [ id  " + i.ToString());
        glycanGraphFile.WriteLine("    LabelGraphics");
        glycanGraphFile.WriteLine("    [ text \"" + label + uniqueGlycansGlycopepList[i].glycanMass + "\" ]");

        if (label != "")
        {
            glycanGraphFile.WriteLine("    graphics");
            glycanGraphFile.WriteLine("    [ fill	\"#4287F5\" ]");
        }
        glycanGraphFile.WriteLine("  ]");

        otherGraph.Add(uniqueGlycansGlycopepList[i], new List<Edge<Glycopep>>());
    }

    for (int i = 0; i < uniqueGlycansGlycopepList.Count; i++)
    {
        for (int j = i + 1; j < uniqueGlycansGlycopepList.Count; j++)
        {
            int biggerIndex, smallerIndex;

            if (uniqueGlycansGlycopepList[i].glycanMass > uniqueGlycansGlycopepList[j].glycanMass)
            {
                biggerIndex = i;
                smallerIndex = j;
            }
            else
            {
                biggerIndex = j;
                smallerIndex = i;
            }

            foreach (string delta in deltas)
            {
                glycanGraphFile.Write(glycans.writeEdge(biggerIndex, smallerIndex, Constants.glycan_masses[delta], delta, uniqueGlycansGlycopepList, otherGraph));
            }
        }
    }

    glycanGraphFile.WriteLine("]");
}

// Preferred way is now based on separating by peptide
for (int k = 0; k < mapPeptoGlycanList.Count; k++)
{
    double keyVal = mapPeptoGlycanList.Keys.ToArray()[k];
    using (var glycanGraphSpecificPeptideFile = new StreamWriter(Path.Combine(spectraDirectory, @"graphs - " + spectraFileName + @"\glycan graph pep - " + keyVal + " - " + spectraFileName + ".gml")))
    {
        int len = mapPeptoGlycanList[keyVal].Count;

        glycanGraphSpecificPeptideFile.WriteLine("Creator \"yFiles\"\nVersion 2.2\ngraph\n[ hierarchic  1\n  directed  1  ");

        for (int i = 0; i < len; i++)
        {
            string label = "";
            if (mapPeptoGlycanList[keyVal][i].comp != null)
            {
                label = mapPeptoGlycanList[keyVal][i].comp.toString() + ": ";
            }
            glycanGraphSpecificPeptideFile.WriteLine("  node");
            glycanGraphSpecificPeptideFile.WriteLine("  [ id  " + i.ToString());
            glycanGraphSpecificPeptideFile.WriteLine("    LabelGraphics");
            glycanGraphSpecificPeptideFile.WriteLine("    [ text \"" + label + mapPeptoGlycanList[keyVal][i].glycanMass + "\nrt: " + mapPeptoGlycanList[keyVal][i].rt + "\" ]");

            if (mapPeptoGlycanList[keyVal][i].comp != null)//(!label.Equals(""))
            {
                glycanGraphSpecificPeptideFile.WriteLine("    graphics");
                glycanGraphSpecificPeptideFile.WriteLine("    [ fill	\"#4287F5\" ]");
            }
            glycanGraphSpecificPeptideFile.WriteLine("  ]");

            //Glycopep key = new Glycopep();
            //key.glycanMass = Convert.ToDouble(mapPeptoGlycanList[keyVal][i].glycanMass);
            //key.peptideMass = keyVal;
            //key.comp = GlycanComposition.toComp(label);

            graph3.Add(mapPeptoGlycanList[keyVal][i], new List<Edge<Glycopep>>());
        }
        
        for (int i = 0; i < len; i++)
        {
            for (int j = i + 1; j < len; j++)
            {
                int biggerIndex, smallerIndex;

                if (mapPeptoGlycanList[keyVal][i].glycanMass > mapPeptoGlycanList[keyVal][j].glycanMass)
                {
                    biggerIndex = i;
                    smallerIndex = j;
                }
                else
                {
                    biggerIndex = j;
                    smallerIndex = i;
                }

                List<Glycopep> list = mapPeptoGlycanList[keyVal];

                foreach (string delta in deltas)
                {
                    glycanGraphSpecificPeptideFile.Write(glycans.writeEdge(biggerIndex, smallerIndex, Constants.glycan_masses[delta], delta, list, graph3));
                }
            }
        }
    }
}

/*
// All-in-one graph is deprecated
List<string> additionalGlycansDeprecated = glycans.FindAdditionalGlycans(glycopepList, graph2, compositions, Path.Combine(spectraDirectory, spectraFileName + " - deprecated - additional glycans - details.csv"));
using (var uniqueAdditionalGlycansFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - deprecated - additional glycans - summary.txt")))
{
    foreach (string g in additionalGlycansDeprecated)
    {
        uniqueAdditionalGlycansFile.WriteLine(g);
    }
}
*/

// Reset - i.e., undo changes made by FindAdditionalGlycans()
foreach (Glycopep glycopep in glycopepList)
{
    if (glycopep.compChanged)
    {
        glycopep.comp = null;
        glycopep.compChanged = false;
    }
}

// Preferred way is now based on separating by peptide
List<string> additionalGlycans = glycans.FindAdditionalGlycans(glycopepList, graph3, compositions, Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - details.csv"));
Console.WriteLine("Additional glycans:");
using (var uniqueAdditionalGlycansFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - summary.txt")))
{
    using (var byonicModificationsFineControlFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - Byonic modifications fine control.txt")))
    {
        using (var byonicColumnFilterFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - Byonic column filter.cft")))
        {
            byonicColumnFilterFile.Write("""{"application":"byos","author":"GlycoNetworks","column_filters":[{"column":"","id":"","operator":"Script","value":"var x = Column.value(\"Glycans\"); return """);
            foreach (string g in additionalGlycans)
            {
                uniqueAdditionalGlycansFile.WriteLine(g);
                byonicModificationsFineControlFile.WriteLine(g + " @ NGlycan | rare1");
                byonicColumnFilterFile.Write(@"x == \""" + g + @"\"" || ");
                Console.WriteLine(g);
            }
            byonicColumnFilterFile.WriteLine("""false;"}],"created":"2025-02-14T15:21:00","name":"Additional glycans filter","table":"Peptides","textFilter":"","type":"Standard","version":1}""");
        }
    }
}

/*
using (var file6 = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName) + " - mostCommonGlycans.csv"))
{
    file6.WriteLine("Glycan Mass, # Occurences, Retention Times of Occurences, Composition, Cluster Size");
    for (int i = 0; i < commonGlycans.Count; i++)
    {
        file6.WriteLine(commonGlycans[i].glycMass + ", " + commonGlycans[i].occurrences + ", " + commonGlycans[i].rtOcc + ", " + commonGlycans[i].comp + ", " + commonGlycans[i].clusterSize);
    }
}
*/

// Alyssa's code starts here
/* FYI: 
 * mapPeptoGlycanList is a Dictionary<double, List<Glycopep>>();
 * var graph3 = new Dictionary<Glycopep, List<Edge>>(); (graph3 is the preferred way of separating graphs based on peptide (aka what is currently used)
 *    for each glycan mass, list of edges connected to it? 
 * edge created btwn two nodes (glycans) IF: ((glycan_mass - Constants.tolerance < diff && diff < glycan_mass + Constants.tolerance))
 *    aka if the difference between the two glycan masses is some glycan mass (within tolerable range) 
 * var mapPeptoGlycanList = new Dictionary<double, List<Glycopep>>();
 */

// VARIABLES & OUTPUT FILES
var DiffValuesFilteredFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - peptide to difference values.txt"));
var initialDiffValuesFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - initial difference value filter for potential blocks.txt"));
var testerFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - tester file.txt"));
var massValueDataFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - mass value data file.csv"));

var mapPepToGlycanGraph = new Dictionary<double, List<GlycoGraph>>();
var mapPepToGraphToDiffValuesSet = new Dictionary<double, Dictionary<Tuple<int, int>, Dictionary<double, List<double>>>>();
var mapPepToGraphsToDiffValues = new Dictionary<double, Dictionary<Tuple<int, int>, List<DiffValueGraph>>>();

foreach (double pep in mapPeptoGlycanList.Keys)
{
    // PRE-PROCESSING
    mapPepToGlycanGraph[pep] = new List<GlycoGraph>();
    mapPepToGraphToDiffValuesSet[pep] = new Dictionary<Tuple<int, int>, Dictionary<double, List<double>>>();
    mapPepToGraphsToDiffValues[pep] = new Dictionary<Tuple<int, int>, List<DiffValueGraph>>();

    int len = mapPeptoGlycanList[pep].Count;

    // initializing each glycan as a small graph with one node. 
    for (int i = 0; i < len; i++)
    {
        GlycoGraph graph = new GlycoGraph();
        graph.numNodes = 1;
        graph.numEdges = 0;
        graph.nodes = new List<Glycopep>();
        graph.nodes = new List<Glycopep>() { mapPeptoGlycanList[pep][i] };
        graph.edges = new List<Edge<Glycopep>>();
        mapPepToGlycanGraph[pep].Add(graph);
    }

    //creating edges + merging as appropriate
    for (int i = 0; i < len; i++)
    {
        for (int j = i+1; j < len; j++)
        {
            Glycopep node1 = mapPeptoGlycanList[pep][i];
            Glycopep node2 = mapPeptoGlycanList[pep][j];
            double diff; // difference between two glycan nodes
            if (node1.glycanMass > node2.glycanMass)
            {
                diff = node1.glycanMass - node2.glycanMass;
            }
            else
            {
                diff = node2.glycanMass - node1.glycanMass;
            }

            foreach (string delta in deltas)
            {
                if ((Constants.glycan_masses[delta] - Constants.tolerance < diff && diff < Constants.glycan_masses[delta] + Constants.tolerance)) // if the difference is approx the value of a known delta, merge the graphs
                {
                    int graph1Index = 0, graph2Index = 0;
                    for (int k=0; k < mapPepToGlycanGraph[pep].Count; k++) // identifying the graphs associated with these glycan nodes
                    {
                        if (mapPepToGlycanGraph[pep][k].nodes.Contains(node1))
                        {
                            graph1Index = k;
                        }
                        if (mapPepToGlycanGraph[pep][k].nodes.Contains(node2))
                        {
                            graph2Index = k;
                        }
                    }

                    if (graph1Index == graph2Index) // if glycan nodes in same graph, add an edge
                    {
                        mapPepToGlycanGraph[pep][graph1Index].addEdge(node1, node2, delta);
                    }
                    else // else, merge the graphs
                    {
                        mapPepToGlycanGraph[pep][graph1Index].mergeGraph(mapPepToGlycanGraph[pep][graph2Index], node1, node2, delta);
                        mapPepToGlycanGraph[pep].RemoveAt(graph2Index);

                    }
                }
            }

        }
    }
    // sort graphs by numNodes
    mapPepToGlycanGraph[pep] = mapPepToGlycanGraph[pep].OrderByDescending(graph => graph.numNodes).ToList();

    // DISCOVERY STEP 
    // brute-force counts the differences between glycan masses in different graphs (for each possible pair), counting how many times that difference appears
    // in other words, diff-value generation 
    var mapGraphsToDiffValues = new Dictionary<Tuple<int, int>, DiffValueBuckets>();
    for (int i=0; i < mapPepToGlycanGraph[pep].Count; i++)
    {
        for (int j=i+1; j < mapPepToGlycanGraph[pep].Count; j++)
        {
            if (mapPepToGlycanGraph[pep][i].numNodes > 2 && mapPepToGlycanGraph[pep][j].numNodes > 2)
            {
                var tuple = Tuple.Create(i, j); // for every tuple (i,j), the diff vlaue = graphJ nodes - graphI nodes 
                mapPepToGraphToDiffValuesSet[pep][tuple] = new Dictionary<double, List<double>>();

                GlycoGraph graphJ = mapPepToGlycanGraph[pep][j]; // graphJ
                GlycoGraph graphI = mapPepToGlycanGraph[pep][i]; // graphI 
                
                // new method of generating + sorting diffValues 
                // note: graph1 is graphJ, and graph2 is graphI 
                List<DiffValue> diffValuesList = new List<DiffValue>();
                foreach (Glycopep nodeJ in graphJ.nodes)
                {
                    foreach (Glycopep nodeI in graphI.nodes)
                    {
                        double diff = nodeJ.glycanMass - nodeI.glycanMass; // graphJ node - graphI node 
                        Tuple<double, Glycopep, Glycopep> nodeTuple = new Tuple<double, Glycopep, Glycopep>(diff, nodeI, nodeJ); 
                        DiffValue diffValue = new DiffValue(diff, graphI, graphJ, nodeI, nodeJ );
                        diffValuesList.Add(diffValue);
                    }
                }

                // ... really this should be in another method or smth.... but lol im lwk lazy? so skullj. 
                Queue<DiffValueBuckets> diffValueBucketsIter = new Queue<DiffValueBuckets>();
                int capacity = 2;
                DiffValueBuckets presentBucket = new DiffValueBuckets(diffValuesList);
                DiffValueBuckets previousBucket = null;
                diffValueBucketsIter.Enqueue(presentBucket);
                Boolean continueIter = true;
                while (continueIter)
                {
                    if (previousBucket != null && (diffValueBucketsIter.Peek().numBuckets <= presentBucket.numBuckets))
                    {
                        continueIter = false;
                    }

                    List<double> medians = presentBucket.returnMedians();
                    previousBucket = presentBucket;
                    presentBucket = new DiffValueBuckets(diffValuesList, medians);
                    diffValueBucketsIter.Enqueue(presentBucket);

                    if (diffValueBucketsIter.Count > capacity)
                    {
                        diffValueBucketsIter.Dequeue();
                    }
                }

                mapGraphsToDiffValues[new Tuple<int, int>(i, j)] = presentBucket;
                
            }
        }
    }

    // sorting diff values into graphs 
    Dictionary<Tuple<int, int>, List<DiffValueGraph>> mapGlycoGraphsToDiffValuesGraph = new Dictionary<Tuple<int, int>, List<DiffValueGraph>>();
    foreach (Tuple<int, int> tuple in mapGraphsToDiffValues.Keys)
    {
        mapGlycoGraphsToDiffValuesGraph[tuple] = new List<DiffValueGraph>();

        // initializing with each diff value as its own subgraph 
        foreach (double repDiffValue in mapGraphsToDiffValues[tuple].buckets.Keys)
        {
            DiffValueGraph graph = new DiffValueGraph();
            graph.numNodes = 1;
            graph.numEdges = 0;

            Tuple<double, Bucket<DiffValue>> node = new Tuple<double, Bucket<DiffValue>>(repDiffValue, mapGraphsToDiffValues[tuple].buckets[repDiffValue]);
            graph.nodes = new List<Tuple<double, Bucket<DiffValue>>>() { node };
            graph.edges = new List<Edge<Tuple<double, Bucket<DiffValue>>>>();
            graph.largestNode = node;
            mapGlycoGraphsToDiffValuesGraph[tuple].Add(graph);
        }

        // merging graphs 
        double[] repDiffValues = mapGraphsToDiffValues[tuple].buckets.Keys.ToArray();
        for (int i=0; i < repDiffValues.Length; i++)
        {
            for (int j=i+1; j < repDiffValues.Length; j++)
            {
                foreach (string delta in deltas)
                {
                    Tuple<double, Bucket<DiffValue>> node1 = new Tuple<double, Bucket<DiffValue>>(repDiffValues[i], mapGraphsToDiffValues[tuple].buckets[repDiffValues[i]]);
                    Tuple<double, Bucket<DiffValue>> node2 = new Tuple<double, Bucket<DiffValue>>(repDiffValues[i], mapGraphsToDiffValues[tuple].buckets[repDiffValues[j]]);
                    double diff;
                    if (repDiffValues[i] > repDiffValues[j]) { diff = repDiffValues[i] - repDiffValues[j]; }
                    else { diff = repDiffValues[j] - repDiffValues[i]; }

                    if ((Constants.glycan_masses[delta] - Constants.tolerance < diff && diff < Constants.glycan_masses[delta] + Constants.tolerance))
                    {
                        int graph1Index = 0, graph2Index = 0;

                        for (int graphIndex = 0; graphIndex < mapGlycoGraphsToDiffValuesGraph[tuple].Count; graphIndex++) // figuring out which graph contains node1 & node2
                        {
                            if (mapGlycoGraphsToDiffValuesGraph[tuple][graphIndex].nodes.Contains(node1)) { graph1Index = graphIndex; }
                            if (mapGlycoGraphsToDiffValuesGraph[tuple][graphIndex].nodes.Contains(node2)) { graph2Index = graphIndex; }
                        }

                        if (graph1Index == graph2Index) // adding edge 
                        {
                            mapGlycoGraphsToDiffValuesGraph[tuple][graph1Index].addEdge(node1, node2, delta);
                        }
                        else // merging two subgraphs 
                        {
                            mapGlycoGraphsToDiffValuesGraph[tuple][graph1Index].mergeGraph(mapGlycoGraphsToDiffValuesGraph[tuple][graph2Index], node1, node2, delta);
                            mapGlycoGraphsToDiffValuesGraph[tuple].RemoveAt(graph2Index);
                        }
                    }
                }
            }
        }
        
        List<DiffValueGraph> sortedDiffValuesGraph = mapGlycoGraphsToDiffValuesGraph[tuple].OrderByDescending(graph => graph.numNodes).ToList();
        // now have graphs of sortedDiffValues!! graphs sorted from most # of nodes to least. 
        mapPepToGraphsToDiffValues[pep][tuple] = sortedDiffValuesGraph;
    }
}

// storing the delta candidates generated from the discovery step -- separated by peptide and by DiffValueGraphs
var sortedDiffValuesGraphs = new List<DiffValueGraph>();
foreach (double pep in mapPepToGraphsToDiffValues.Keys)
{
    DiffValuesFilteredFile.WriteLine(pep + ": ");
    foreach (Tuple<int, int> tuple in mapPepToGraphsToDiffValues[pep].Keys)
    {
        DiffValuesFilteredFile.Write("\t" + "i = " + tuple.Item1 + ", j = " + tuple.Item2 + ": [");
        foreach (DiffValueGraph graph in mapPepToGraphsToDiffValues[pep][tuple])
        {
            DiffValuesFilteredFile.Write("(numNodes = " + graph.numNodes + ", largestNode: ("+ graph.largestNode.Item1 + ", " + graph.largestNode.Item2.items.Count + ")), "); // effectively showing the size + the largestNode of each DiffValueGraph (under those specific i/j) 
            sortedDiffValuesGraphs.Add(graph);
        }
    }
}

sortedDiffValuesGraphs = sortedDiffValuesGraphs.OrderByDescending(item => item.largestNode.Item2.items.Count).ToList();

// storing original delta candidates generated from the discovery step
// also stores "repeated" candidates -- candidates that differ by a building block, such as a Hex or HexNAc 
initialDiffValuesFile.WriteLine("Potential building blocks: ");
for (int i = 0; i<20 && i<sortedDiffValuesGraphs.Count; i++)
{
    List<Tuple<double, Bucket<DiffValue>>> sortedGraphNodes = sortedDiffValuesGraphs[i].nodes.OrderByDescending(tuple => tuple.Item2.items.Count).ToList();
    initialDiffValuesFile.WriteLine("Block Mass: " + sortedGraphNodes[0].Item1 + ", Number of Appearances: " + sortedGraphNodes[0].Item2.items.Count);
    for (int j=1; j<10 && j < sortedGraphNodes.Count; j++)
    {
        initialDiffValuesFile.WriteLine("\tUnderling - " + "Block Mass: " + sortedGraphNodes[j].Item1 + ", Number of Appearances: " + sortedGraphNodes[j].Item2.items.Count);
    }
}
initialDiffValuesFile.WriteLine();

// from all of the potential delta candidates generated across differing peptides, select the topmost candidates 
// also take a few "repeated" candidates (i.e. off by a Hex or HexNAc)
// this list of potential delta candidates will then go through the confirmation step (described later) 
var potentialBlocks = new List<Block>();
for (int i = 0; i<100 && i<sortedDiffValuesGraphs.Count; i++)
{
    for (int j = 0; j<5 && j < sortedDiffValuesGraphs[i].numNodes; j++)
    {
        Block block = new Block();

        sortedDiffValuesGraphs[i].nodes = sortedDiffValuesGraphs[i].nodes.OrderByDescending(tuple => tuple.Item2.items.Count).ToList();
        block.diffValueGraph = sortedDiffValuesGraphs[i];
        double mass = sortedDiffValuesGraphs[i].nodes[j].Item1;
        block.diffValueHits = sortedDiffValuesGraphs[i].nodes[j].Item2.items.Count;

        HashSet<Glycopep> graphIRelevantNodes = new HashSet<Glycopep>();
        HashSet<Glycopep> graphJRelevantNodes = new HashSet<Glycopep>();
        foreach (DiffValue diffValue in sortedDiffValuesGraphs[i].nodes[j].Item2.items)
        {
            graphIRelevantNodes.Add(diffValue.nodeGraphI);
            graphJRelevantNodes.Add(diffValue.nodeGraphJ);
        }

        if (mass > 0 )
        {
            block.mass = mass; // means graphJ nodes - graphI nodes is positive
            block.smallerGlycographRelevantNodes = graphIRelevantNodes;
            block.largerGlycographRelevantNodes = graphJRelevantNodes;
        }
        else
        {
            block.mass = -mass;
            block.smallerGlycographRelevantNodes = graphJRelevantNodes;
            block.largerGlycographRelevantNodes = graphIRelevantNodes;
        }

        block.deltaHits = 0;
        block.controlDeltaPeaks = new List<Tuple<MS2, Point>>();
        block.controlDeltaBlockPeaks = new List<Tuple<MS2, Point>>();
        block.deltaPeaks = new List<Tuple<MS2, Point>>();
        block.deltaBlockPeaks = new List<Tuple<MS2, Point>>();

        block.controlDeltaConfidence = 0;
        block.controlDeltaBlockConfidence = 0;
        block.deltaConfidence = 0;
        block.deltaBlockConfidence = 0;
        block.confidence = 0;

        block.realMassDifferences = new List<realMassDifferences>();
        block.realMassValues = new List<double>();
        block.medianRealMassValue = 0;

        potentialBlocks.Add(block);
    }
}

testerFile.Write("deltas: [");
foreach (string delta in deltas)
{
    testerFile.Write(delta + ", ");
}
testerFile.WriteLine("]");

// CONFIRMATION STEP 
// examining ms/ms fragmentation spectra from the glycograph pair that generates the according diffValueGraph & potential block 
// from these ms/ms spectra, finding the (building block mass) peaks and the (building block mass + potential delta candidate mass) peaks 
// these confirmation peaks are used to produce a more precise delta candidate mass 
// the intensity of these peaks is used to calculate the quality score of the hypothetical delta candidate

var glycanIonSubset = new Dictionary<string, double>();

foreach (string delta in deltas)
{
    foreach (string key in Constants.glycanSpecificIons.Keys)
    {
        if (key == delta)
        {
            glycanIonSubset.Add(key, Constants.glycanSpecificIons[key]);
        }
    }
}
if (glycanIonSubset.ContainsKey("Hex-HexNAc")) glycanIonSubset.Add("Hex(2)-HexNAc", glycanIonSubset["Hex-HexNAc"] + Constants.glycan_masses["Hex"]);

double HexHexNAc = Constants.glycan_masses["Hex-HexNAc"] + Constants.adduct;
foreach (Block block in potentialBlocks)
{
    block.numSpectra = 0;
    block.spectra = new List<MS2>();
    block.indices = new List<int>();
    block.RTToMass = new Dictionary<MS2, List<double>>();
    foreach (MS2 spectrum in spectra)
    {
        
        // DEPRECATED CODE 
        foreach (Glycopep pep in block.smallerGlycographRelevantNodes)
        {
            if (spectrum.PrecursorMz == pep.precursorMz) // effectively going thru spectra that generate the glycographs 
            {
                // expecting Hex-HexNAc (366) peak in spectra
                Point HexHexNAcPeak = null;
                foreach (Point datapoint in spectrum.getData())
                {
                    if (HexHexNAc - Constants.tolerance < datapoint.mass && datapoint.mass < HexHexNAc + Constants.tolerance)
                    {
                        HexHexNAcPeak = datapoint;
                    }
                }
                if (HexHexNAcPeak == null)
                {
                    continue;
                }

                var spectrumBlocks = new Dictionary<double, Boolean>();
                var spectrumBlocksTrueValues = new Dictionary<double, double>();
                foreach (string delta in glycanIonSubset.Keys)
                {
                    spectrumBlocks[glycanIonSubset[delta]] = false;
                    if (delta == "NeuAc" || delta == "NeuGc")
                    {
                        spectrumBlocks[glycanIonSubset[delta] - Constants.water] = false;
                        // spectrumBlocks[glycanIonSubset[delta] - 2 * Constants.water] = false;
                    }
                }

                foreach (Point datapoint in spectrum.getData())
                {
                    foreach (double delta in spectrumBlocks.Keys)
                    {
                        if (delta - Constants.tolerance < datapoint.mass && datapoint.mass < delta + Constants.tolerance)
                        {
                            spectrumBlocks[delta] = true;
                            block.controlDeltaPeaks.Add(new Tuple<MS2, Point>(spectrum, datapoint));
                            spectrumBlocksTrueValues[delta] = datapoint.mass;

                            // calculating confidence (standardized relative to Hex-HexNAc, mass 366) 
                            double variation = 1 - datapoint.intensity / HexHexNAcPeak.intensity;
                            block.controlDeltaConfidence += 1 - Constants.weight * variation;
                        }
                    }
                }

                foreach (Point datapoint in spectrum.getData())
                {
                    foreach (string delta in glycanIonSubset.Keys)
                    {
                        if (delta != "NeuAc" && delta != "NeuGc")
                        {
                            double deltaMass = glycanIonSubset[delta];
                            if (spectrumBlocks[deltaMass]) // no use looking for delta (known building block) mass + hypothetical building block mass if the relevant delta peak doesn't exist in the spectra 
                            {
                                double trueDeltaMass = spectrumBlocksTrueValues[glycanIonSubset[delta]];
                                Boolean deltaBlockNull = false; // checking to see if delta block.mass is also equivalent/ballpark another delta, and hence should not be counted 
                                foreach (double delta2 in spectrumBlocks.Keys)
                                {
                                    if (trueDeltaMass + block.mass - Constants.tolerance < delta2 && delta2 < trueDeltaMass + block.mass + Constants.tolerance) { deltaBlockNull = true; }
                                }

                                if (deltaBlockNull == false && (trueDeltaMass + block.mass - Constants.tolerance < datapoint.mass && datapoint.mass < trueDeltaMass + block.mass + Constants.tolerance)) // additionally checking if datapoint mass is within tolerance range of deltaMass + blockMass 
                                {
                                    block.controlDeltaBlockPeaks.Add(new Tuple<MS2, Point>(spectrum, datapoint));

                                    // calculating confidence (standardized relative to Hex-HexNAc, mass 366) 
                                    double variation = 1 - datapoint.intensity / HexHexNAcPeak.intensity;
                                    block.controlDeltaBlockConfidence += 1 - Constants.weight * variation;
                                }
                            }
                        }
                        else
                        {
                            double deltaMass = glycanIonSubset[delta];
                            if (spectrumBlocks[deltaMass] && spectrumBlocks[deltaMass - Constants.water])
                            {
                                double trueDeltaMass = spectrumBlocksTrueValues[glycanIonSubset[delta]];
                                Boolean deltaBlockNull = false; // checking to see if delta block.mass is also equivalent/ballpark another delta, and hence should not be counted 
                                foreach (double delta2 in spectrumBlocks.Keys)
                                {
                                    if (trueDeltaMass + block.mass - Constants.tolerance < delta2 && delta2 < trueDeltaMass + block.mass + Constants.tolerance) { deltaBlockNull = true; }
                                }

                                if (deltaBlockNull == false && (trueDeltaMass + block.mass - Constants.tolerance < datapoint.mass && datapoint.mass < trueDeltaMass + block.mass + Constants.tolerance)) // additionally checking if datapoint mass is within tolerance range of deltaMass + blockMass 
                                {
                                    block.controlDeltaBlockPeaks.Add(new Tuple<MS2, Point>(spectrum, datapoint));

                                    // calculating confidence (standardized relative to Hex-HexNAc, mass 366) 
                                    double variation = 1 - datapoint.intensity / HexHexNAcPeak.intensity;
                                    block.controlDeltaBlockConfidence += 1 - Constants.weight * variation;
                                }   
                            }
                        }
                    }
                }
            }
        }
        // END OF DEPRECATED CODE

        bool addSpectraCount = false;
        var masses = new List<double>();

        // examining the spectra associated with the glycan network nodes with the +delta mass 
        foreach (Glycopep pep in block.largerGlycographRelevantNodes)
        {
            if (spectrum.PrecursorMz == pep.precursorMz)
            {
                // expecting Hex-HexNAc (366) peak in spectra
                // using this peak for quality score (sum of peak intensities relative to Hex-HexNAc peak)
                Point HexHexNAcPeak = null;
                foreach (Point datapoint in spectrum.getData())
                {
                    if (HexHexNAc - Constants.tolerance / 5 < datapoint.mass && datapoint.mass < HexHexNAc + Constants.tolerance / 5)
                    {
                        HexHexNAcPeak = datapoint;
                    }
                }
                if (HexHexNAcPeak == null)
                {
                    continue;
                }

                // recording oxionium ion peaks/building block peaks in the spectra 
                var spectrumBlocks = new Dictionary<double, Boolean>();
                var spectrumBlocksTrueValues = new Dictionary<double, double>();
                foreach (string delta in glycanIonSubset.Keys)
                {
                    spectrumBlocks[glycanIonSubset[delta]] = false;
                    if (delta == "NeuAc" || delta == "NeuGc")
                    {
                        spectrumBlocks[glycanIonSubset[delta] - Constants.water] = false;
                    }
                }

                foreach (Point datapoint in spectrum.getData()) 
                {
                    foreach (double delta in spectrumBlocks.Keys)
                    {
                        if (delta - Constants.tolerance / 5 < datapoint.mass && datapoint.mass < delta + Constants.tolerance / 5)
                        {
                            spectrumBlocks[delta] = true;
                            block.deltaPeaks.Add(new Tuple<MS2, Point>(spectrum, datapoint));
                            spectrumBlocksTrueValues[delta] = datapoint.mass;

                            // calculating confidence (standardized relative to Hex-HexNAc, mass 366) 
                            double variation = 1 - datapoint.intensity / HexHexNAcPeak.intensity;
                            block.deltaConfidence += 1 - Constants.weight * variation;
                        }
                    }
                }

                // investigating datapoints with mass (building block mass + potential candidate mass)
                foreach (Point datapoint in spectrum.getData())
                {
                    foreach (string delta in glycanIonSubset.Keys)
                    {
                        double deltaMass = glycanIonSubset[delta];
                        Boolean sialicAcid = (delta == "NeuAc" || delta == "NeuGc");
                        if (spectrumBlocks[deltaMass] && (!sialicAcid || spectrumBlocks[deltaMass - Constants.water])) // no use looking for delta (known building block) mass + hypothetical building block mass if the relevant delta peak doesn't exist in the spectra 
                        {
                            double trueDeltaMass = spectrumBlocksTrueValues[glycanIonSubset[delta]];
                            Boolean deltaBlockNull = false; // checking to see if delta block.mass is also equivalent/ballpark another delta, and hence should not be counted 
                            foreach (double delta2 in spectrumBlocks.Keys)
                            {
                                if (trueDeltaMass + block.mass - Constants.tolerance / 5 < delta2 && delta2 < trueDeltaMass + block.mass + Constants.tolerance / 5) { deltaBlockNull = true; }
                            }

                            if (deltaBlockNull == false && (trueDeltaMass + block.mass - Constants.tolerance / 2 < datapoint.mass && datapoint.mass < trueDeltaMass + block.mass + Constants.tolerance / 2)) // additionally checking if datapoint mass is within tolerance range of deltaMass + blockMass 
                            {
                                block.deltaHits++;
                                block.deltaBlockPeaks.Add(new Tuple<MS2, Point>(spectrum, datapoint));
                                block.realMassValues.Add(datapoint.mass - trueDeltaMass);
                                realMassDifferences difference = new realMassDifferences(datapoint.mass - trueDeltaMass, trueDeltaMass, datapoint.mass, spectrum);
                                block.realMassDifferences.Add(difference);

                                // calculating quality score aka confidence (standardized relative to Hex-HexNAc, mass 366) 
                                double variation = 1 - datapoint.intensity / HexHexNAcPeak.intensity;
                                block.deltaBlockConfidence += 1 - Constants.weight * variation;
                                    
                                addSpectraCount = true;
                                masses.Add(datapoint.mass - deltaMass);
                            }
                        }
                    }        
                }
            }
        }

        if (addSpectraCount)
        {
            block.numSpectra++;
            block.spectra.Add(spectrum);
            block.indices.Add(block.realMassValues.Count - 1);
            block.RTToMass[spectrum] = masses;
        }
    }
}

double confidenceWeight = 0; 
foreach (Block block in potentialBlocks)
{
    block.confidence = (1-confidenceWeight) * block.deltaBlockConfidence - (confidenceWeight) * block.controlDeltaBlockConfidence;
    // originally, confidence took into account both deltaBlockConfidence and controlDeltaBlockConfidence (see deprecated code) 
    // now, quality score/confidence is just the deltaBlock confidence -- the sum of the confirmation peak intensities standardized to the Hex-HexNAc peak
    // this is due to suspected leakage because precursor m/z isolation window widths are not zero, leading to unusually high controDeltaBlockConfidence for some datasets
    //     (for instance, with the old metric, the Hi-5 dataset didn't show phosphorylcholine as a candidate) 
}

// identifying median value for the confirmation peaks of each potential delta
foreach (Block block in potentialBlocks)
{
    block.realMassValues = block.realMassValues.Order().ToList();
    int medianIndex = block.realMassValues.Count / 2;
    if (block.realMassValues.Count > 0)
    {
        if (block.realMassValues.Count % 2 == 0) { block.medianRealMassValue = (block.realMassValues[medianIndex-1] + block.realMassValues[medianIndex])/2; }
        else { block.medianRealMassValue = block.realMassValues[medianIndex]; }
    }

    // delta mass determined to be the median of realMassValues
    if (block.medianRealMassValue != 0)
    {
        block.mass = block.medianRealMassValue;
    }
}
    
// Storing data on each potential delta candidate before merging similar candidates 
// Note: original potential deltas were selected from candidates generated from separate peptides (GlycoNetworks separates glycan networks by peptide)
//       hence, effectively the same delta candidate is bound to arise from different peptides as well 
foreach (Block block in potentialBlocks)
{
    testerFile.WriteLine(block.mass + ": ");
    //testerFile.WriteLine("'real' median mass value: " + block.medianRealMassValue);
    testerFile.WriteLine("diffValueHits = " + block.diffValueHits);
    testerFile.WriteLine("deltaHits = " + block.deltaHits);

    block.controlDeltaPeaks = block.controlDeltaPeaks.OrderByDescending(tuple => tuple.Item2.intensity).ToList();
    block.controlDeltaBlockPeaks = block.controlDeltaBlockPeaks.OrderByDescending(tuple => tuple.Item2.intensity).ToList();
    block.deltaPeaks = block.deltaPeaks.OrderByDescending(tuple => tuple.Item2.intensity).ToList();
    block.deltaBlockPeaks = block.deltaBlockPeaks.OrderByDescending((tuple) => tuple.Item2.intensity).ToList();

    testerFile.Write("controlDeltaPeaks = [");
    foreach (Tuple<MS2, Point> peak in block.controlDeltaPeaks)
    {
        testerFile.Write("(precursorMz = " + peak.Item1.PrecursorMz + ", mass = " + peak.Item2.mass + ", intensity = " + peak.Item2.intensity + "), ");
    }
    testerFile.WriteLine("]");

    testerFile.Write("controlDeltaBlockPeaks = [");
    foreach (Tuple<MS2, Point> peak in block.controlDeltaBlockPeaks)
    {
        testerFile.Write("(precursorMz = " + peak.Item1.PrecursorMz + ", mass = " + peak.Item2.mass + ", intensity = " + peak.Item2.intensity + "), ");
    }
    testerFile.WriteLine("]");

    testerFile.Write("deltaPeaks = [");
    foreach (Tuple<MS2, Point> peak in block.deltaPeaks)
    {
        testerFile.Write("(precursorMz = " + peak.Item1.PrecursorMz + ", mass = " + peak.Item2.mass + ", intensity = " + peak.Item2.intensity + "), ");
    }
    testerFile.WriteLine("]");

    testerFile.Write("deltaBlockPeaks = [");
    foreach (Tuple<MS2, Point> peak in block.deltaBlockPeaks)
    {
        testerFile.Write("(precursorMz = " + peak.Item1.PrecursorMz + ", mass = " + peak.Item2.mass + ", intensity = " + peak.Item2.intensity + "), ");
    }
    testerFile.WriteLine("]");

    testerFile.WriteLine("Confidences: ");
    testerFile.WriteLine("\tcontrolDeltaConfidence = " + block.controlDeltaConfidence);
    testerFile.WriteLine("\tcontrolDeltaBlockConfidence = " + block.controlDeltaBlockConfidence);
    testerFile.WriteLine("\tdeltaConfidence = " + block.deltaConfidence);
    testerFile.WriteLine("\tdeltaBlockConfidence = " + block.deltaBlockConfidence);
    testerFile.WriteLine("\tblock confidence (aka deltaBlockConfidence - controlDeltaBlockConfidence) = " + block.confidence);

    testerFile.WriteLine();
}

potentialBlocks = potentialBlocks.OrderByDescending(block => block.confidence).ToList();

// merging similar delta candidates into one BlockBucket via the BlockBucketSorter 
var bucketIterations = new Queue<BlockBuckets>();
BlockBuckets currentBucket = new BlockBuckets(potentialBlocks);
BlockBuckets lastBucket = null;
bucketIterations.Enqueue(currentBucket);
Boolean continueOn = true;
while (continueOn)
{
    List<double> medians = currentBucket.returnMedians();
    lastBucket = currentBucket;
    currentBucket = new BlockBuckets(potentialBlocks, medians);
    bucketIterations.Enqueue(currentBucket);
    bucketIterations.Dequeue();

    if (lastBucket.numBuckets <= currentBucket.numBuckets)
    {
        continueOn = false;
    }
}

// merging real mass values for each BlockBucket
foreach (BlockBucket bucket in currentBucket.blockBuckets.Values)
{
    foreach (Block block in bucket.items)
    {
        foreach (double value in block.realMassValues)
        {
            bucket.mergedRealMassValues.Add(value);
        }
        foreach (realMassDifferences value in block.realMassDifferences)
        {
            bucket.mergedRealMassDifferences.Add(value);
        }
        bucket.mergedConfidence += block.confidence;
        bucket.numSpectraTotal += block.numSpectra;
    }
}

// removing outliers via statistics 
foreach (BlockBucket bucket in currentBucket.blockBuckets.Values)
{
    if (bucket.mergedConfidence > 0 && bucket.mergedRealMassValues.Count > 1)
    {
        bucket.medianValue = MathNet.Numerics.Statistics.Statistics.Median(bucket.mergedRealMassValues);
        double IQR = MathNet.Numerics.Statistics.Statistics.InterquartileRange(bucket.mergedRealMassValues);
        double lowerBound;
        double upperBound;
        if (IQR > 0)
        {
            lowerBound = MathNet.Numerics.Statistics.Statistics.LowerQuartile(bucket.mergedRealMassValues) - 1.5 * IQR;
            upperBound = MathNet.Numerics.Statistics.Statistics.UpperQuartile(bucket.mergedRealMassValues) + 1.5 * IQR;
        }
        else
        {
            lowerBound = MathNet.Numerics.Statistics.Statistics.Minimum(bucket.mergedRealMassValues);
            upperBound = MathNet.Numerics.Statistics.Statistics.Maximum(bucket.mergedRealMassValues);
        }

        // revising # of real mass values according to bounds 
        List<double> newRealMassValues = new List<double>();
        List<realMassDifferences> newRealMassDifferences = new List<realMassDifferences>();
        foreach (double value in bucket.mergedRealMassValues)
        {
            if (lowerBound <= value && value <= upperBound) { newRealMassValues.Add(value);  }
        }
        foreach (realMassDifferences value in bucket.mergedRealMassDifferences)
        {
            if (lowerBound <= value.diff && value.diff <= upperBound) { newRealMassDifferences.Add(value); }
        }
        bucket.mergedRealMassValues = newRealMassValues;
        bucket.mergedRealMassDifferences = newRealMassDifferences;

    }
}

// counting number of samples (aka oxoniumIonMasses)
foreach (BlockBucket bucket in currentBucket.blockBuckets.Values)
{
    if (bucket.mergedConfidence > 0 && bucket.mergedRealMassValues.Count > 1)
    {
        foreach (realMassDifferences diff in bucket.mergedRealMassDifferences)
        {
            int oxoniumIonMass = Convert.ToInt32(diff.smallerValue);
            bucket.oxoniumIonMasses.Add(oxoniumIonMass);
        }
    }
}

// confidence interval for building block masses 
foreach (BlockBucket bucket in currentBucket.blockBuckets.Values)
{
    if (bucket.mergedConfidence > 0 && bucket.mergedRealMassValues.Count > 1)
    {
        bucket.statistics = new DescriptiveStatistics(bucket.mergedRealMassValues);
        //double criticalValue = MathNet.Numerics.Distributions.StudentT.InvCDF(0, 1, bucket.mergedRealMassValues.Count - 1, 0.975);
        int df = bucket.oxoniumIonMasses.Count - 1;
        if (df < 2) df = 2;
        double criticalValue = MathNet.Numerics.Distributions.StudentT.InvCDF(0, 1, df, 0.975);
        //Console.Write(criticalValue + "-tester:" + MathNet.Numerics.Distributions.StudentT.InvCDF(0, 1, 70, 0.995) + ", ");
        bucket.marginOfError = criticalValue * bucket.statistics.StandardDeviation / System.Math.Sqrt(bucket.oxoniumIonMasses.Count);
        //bucket.marginOfError = criticalValue * bucket.statistics.StandardDeviation;
        bucket.marginOfError2 = criticalValue * bucket.statistics.StandardDeviation;
        bucket.medianValue = MathNet.Numerics.Statistics.Statistics.Median(bucket.mergedRealMassValues);
    }
}

// putting potential deltas into a network once more to remove "repeat" candidates (aka, those that are simply a Hex or HexNAc off) 
// initializing mergedBlocks from each delta candidate (BlockBucket from BlockBucketSorter)
var mergedBlocks = new List<MergedBlock>();
foreach (double key in currentBucket.blockBuckets.Keys)
{
    MergedBlock block = new MergedBlock();
    block.mass = currentBucket.blockBuckets[key].medianValue;
    block.realMassValues = currentBucket.blockBuckets[key].mergedRealMassValues;
    block.confidence = currentBucket.blockBuckets[key].mergedConfidence;

    block.blockBucket = currentBucket.blockBuckets[key];
    block.keyValuePair = new Tuple<double, BlockBucket>(key, block.blockBucket);

    mergedBlocks.Add(block);
}

List<MergedBlockGraph> blockGraph = new List<MergedBlockGraph>();

// initializing subgraphs 
foreach (MergedBlock block in mergedBlocks)
{
    if (block.mass > 0)
    {
        MergedBlockGraph graph = new MergedBlockGraph();
        graph.numNodes = 1;
        graph.numEdges = 0;

        graph.nodes = new List<MergedBlock>() { block };
        graph.edges = new List<Edge<MergedBlock>>();
        graph.largestNode = block;
        blockGraph.Add(graph);
    }
}

// merging subgraphs
for (int i=0; i<mergedBlocks.Count; i++)
{
    for (int j=i+1; j<mergedBlocks.Count; j++)
    {
        foreach (string delta in deltas)
        {
            MergedBlock node1 = mergedBlocks[i];
            MergedBlock node2 = mergedBlocks[j];

            double diff;
            if (node1.mass > node2.mass) { diff = node1.mass - node2.mass; }
            else { diff = node2.mass - node1.mass;  }

            if ((Constants.glycan_masses[delta] - Constants.tolerance < diff && diff < Constants.glycan_masses[delta] + Constants.tolerance))
            {
                int graph1Index = 0, graph2Index = 0;

                for (int graphIndex = 0; graphIndex < blockGraph.Count; graphIndex++) // figuring out which graph contains node1 & node2
                {
                    if (blockGraph[graphIndex].nodes.Contains(node1)) { graph1Index = graphIndex; }
                    if (blockGraph[graphIndex].nodes.Contains(node2)) { graph2Index = graphIndex; }
                }

                if (graph1Index == graph2Index) // adding edge 
                {
                    blockGraph[graph1Index].addEdge(node1, node2, delta);
                }
                else // merging two subgraphs 
                {
                    blockGraph[graph1Index].mergeGraph(blockGraph[graph2Index], node1, node2, delta);
                    blockGraph.RemoveAt(graph2Index);
                }
            }
        }
    }
}

// sorting blockGraph by a graph's largestNode's quality score
// for each graph, sorting nodes by quality score as well
blockGraph = blockGraph.OrderByDescending(block => block.largestNode.confidence).ToList();
foreach (MergedBlockGraph graph in blockGraph)
{
    graph.nodes = graph.nodes.OrderByDescending(node => node.confidence).ToList();
}

// OUTPUT(S)
// output to command line 
Console.WriteLine();
Console.WriteLine("Candidate delta masses: ");
foreach (MergedBlockGraph graph in blockGraph)
{
    bool suspectedIsotopeArtifact = false;
    for (int i = 0; i < graph.numNodes && i < 2; i++)
    {
        if (Math.Abs(graph.nodes[i].mass - 1.003355) < 0.001) suspectedIsotopeArtifact = true;
    }
    if (suspectedIsotopeArtifact) continue;

    MergedBlock block = graph.largestNode;
    int padding = 3 - (block.mass.ToString("F0").Length);
    string pad = "";
    if (padding > 0) pad = new string(' ', padding);

    if (Constants.verbosityLevel <= 1)
    {
        if (block.confidence < 1 || block.confidence < blockGraph[0].largestNode.confidence / 5) break; // only outputting the top candidates (quality score > 1 and > 1/5 of the top candidate's quality score)
    }

    if (Constants.verbosityLevel >= 0)
    {
        Console.Write("Mass: " + block.mass.ToString("F4") + pad + "\tCI: (" + (block.mass - block.blockBucket.marginOfError).ToString("F4") + ", " + (block.mass + block.blockBucket.marginOfError).ToString("F4") + ")" + pad + pad);
    }

    if (Constants.verbosityLevel >= 1)
    {
        Console.Write("\tQuality score: " + block.confidence.ToString("F1")); // sum of confirmation peak intensities (relative to Hex-HexNAc peak) 
        if (Constants.verbosityLevel >= 2) Console.Write("\tNumber of confirmation peaks: " + block.realMassValues.Count);
    }

    Console.WriteLine();

    if (Constants.verbosityLevel >= 2)
    {
        for (int i = 1; i < graph.numNodes; i++)
        {
            MergedBlock underling = graph.nodes[i];
            Console.Write("\tMass: " + underling.mass.ToString("F4") + pad + "\tCI: (" + (underling.mass - underling.blockBucket.marginOfError).ToString("F4") + ", " + (underling.mass + underling.blockBucket.marginOfError).ToString("F4") + ")" + pad + pad);
            Console.Write("\tQuality score: " + underling.confidence.ToString("F4")); // sum of confirmation peak intensities (relative to Hex-HexNAc peak) 
            Console.Write("\tNumber of confirmation peaks: " + underling.realMassValues.Count);
            Console.WriteLine();
        }
    }
}
 
// gathering metadata on the topmost delta candidate & storing it into massValueDataFile
massValueDataFile.WriteLine("RetentionTime,PrecursorMz,OxoniumMass,IntegerOxoniumMass,OxoniumMassPlusDelta,IntegerOxoniumMassPlusDelta,DeltaMass");
if (blockGraph.Count > 0)
{
    MergedBlock topBlock = blockGraph[0].largestNode;
    foreach (realMassDifferences value in topBlock.blockBucket.mergedRealMassDifferences)
    {
        massValueDataFile.WriteLine(value.spectra.Rt + "," + value.spectra.PrecursorMz + "," + value.smallerValue + "," + value.smallerValue.ToString("F0", System.Globalization.CultureInfo.InvariantCulture) + "," + value.largerValue + "," + value.largerValue.ToString("F0", System.Globalization.CultureInfo.InvariantCulture) + "," + value.diff);
    }
}

DiffValuesFilteredFile.Close();
initialDiffValuesFile.Close();
testerFile.Close();
massValueDataFile.Close();
