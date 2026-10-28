// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>
// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal class CompareGlycan
    {
        public List<Glycopep> visited = new List<Glycopep>();

        public double calcMass(GlycanComposition composition)
        {
            double totalMass = 0;
            for (int i = 0; i < composition.Keys.Count; i++)
            {
                String key = composition.Keys.ToList()[i];
                int number = composition[key];
                double mass = Constants.glycan_masses[key];
                totalMass += number * mass;
            }
            
            return totalMass;
        }

        public string writeEdge(int biggerIndex, int smallerIndex, double glycan_mass, string labelname, List<Glycopep> glycopeptides, Dictionary<Glycopep, List<Edge<Glycopep>>> graph)
        {   
            string edge_code = "";
            double bigger = glycopeptides[biggerIndex].glycanMass;
            double smaller = glycopeptides[smallerIndex].glycanMass;
            double diff = bigger - smaller;
            

            if ((glycan_mass - Constants.tolerance < diff && diff < glycan_mass + Constants.tolerance))
            {
                edge_code = "  edge\n  [ source  " + smallerIndex.ToString() + "\n    target  " + biggerIndex.ToString() + "\n    LabelGraphics" + 
                    "\n    [ text \"" + labelname + "\" ]\n  ]\n";

                Glycopep glycopep1 = glycopeptides[biggerIndex];
                Glycopep glycopep2 = glycopeptides[smallerIndex];
                
                Edge<Glycopep> edge1 = new Edge<Glycopep>();
                edge1.source = glycopep1;
                edge1.target = glycopep2;
                edge1.glycDiff = labelname;
                edge1.increasing = -1;

                Edge<Glycopep> edge2 = new Edge<Glycopep>();
                edge2.source = glycopep2;
                edge2.target = glycopep1;
                edge2.glycDiff = labelname;
                edge2.increasing = 1;

                graph[glycopep1].Add(edge1);
                graph[glycopep2].Add(edge2);
            }

            return edge_code;
        }

        public void dfs(Dictionary<Glycopep, List<Edge<Glycopep>>> graph, Glycopep node, double tolerance, List<Glycopep> visited)
        {
            visited.Add(node);
            //Console.WriteLine(node.glycanMass);
            foreach (Edge<Glycopep> neighbor in graph[node])
            {
                if (!visited.Contains(neighbor.target))
                {
                    //Console.WriteLine("Hello World");
                    dfs(graph, neighbor.target, tolerance, visited);
                }
            }

        }

        public List<string> FindAdditionalGlycans(List<Glycopep> glycopepList, Dictionary<Glycopep, List<Edge<Glycopep>>> graph, List<GlycanComposition> compositions, string additionalGlycansDetailsFilename)
        {
            var additionalGlycans = new List<GlycanComposition>();
            using (var additionalGlycansDetailsFile = new StreamWriter(additionalGlycansDetailsFilename))
            {
                additionalGlycansDetailsFile.WriteLine(@"Glycan composition,Glycan mass,Cluster size,Peptide mass,Retention time");
                bool continueIterating = true;
                while (continueIterating)
                {
                    continueIterating = false;
                    for (int i = 0; i < glycopepList.Count; i++)
                    {
                        if (!graph.ContainsKey(glycopepList[i])) continue;
                        visited.Clear();
                        dfs(graph, glycopepList[i], Constants.tolerance, visited);
                        if (glycopepList[i].comp != null && visited.Count >= 3)
                        {
                            foreach (Edge<Glycopep> edge in graph[glycopepList[i]])
                            {
                                if (edge.target.comp == null)
                                {
                                    continueIterating = true;
                                    edge.target.comp = GlycanComposition.toComp(edge.source.comp.toString());
                                    edge.target.compChanged = true;
                                    if (edge.glycDiff == "Hex-HexNAc")
                                    {
                                        // Hex-HexNAc disaccharide delta special case
                                        if (!edge.target.comp.ContainsKey("Hex")) edge.target.comp.Add("Hex", 0);
                                        if (!edge.target.comp.ContainsKey("HexNAc")) edge.target.comp.Add("HexNAc", 0);
                                        edge.target.comp["Hex"] += edge.increasing;
                                        edge.target.comp["HexNAc"] += edge.increasing;
                                    }
                                    else
                                    {
                                        // Monosaccharide delta
                                        if (!edge.target.comp.ContainsKey(edge.glycDiff)) edge.target.comp.Add(edge.glycDiff, 0);
                                        edge.target.comp[edge.glycDiff] += edge.increasing;
                                    }
                                    if (edge.target.comp.isValid() && !GlycanComposition.ContainsMassMatch(compositions, edge.target.comp) && !GlycanComposition.ContainsMassMatch(additionalGlycans, edge.target.comp))
                                    {
                                        additionalGlycansDetailsFile.WriteLine(edge.target.comp.toString() + "," + edge.target.glycanMass + "," + visited.Count + "," + edge.target.peptideMass + "," + edge.target.rt);
                                        additionalGlycans.Add(edge.target.comp);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            var additionalGlycansAsString = new List<string>();
            foreach (var additionalGlycan in additionalGlycans)
            {
                additionalGlycansAsString.Add(additionalGlycan.toString());
            }
            additionalGlycansAsString.Sort();

            return additionalGlycansAsString;
        }
    }   
}
