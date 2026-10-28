// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>
// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    public class MS2
    {
        List<Point> Data = new List<Point>();
        List<Point> sortedData = new List<Point>();
        int z;

        Point mz;
        Point plusHexnac;
        Point plus2Hexnac;
        Point plusHex2Hexnac;

        public int Charge { get; set; }
        public double Rt { get; set; }
        public double PrecursorMz { get; set; }

        public int findLen()
        {
            return Data.Count;
        }

        public List<Point> getData()
        {
            return Data;
        }

        public double corrIntensity(double mass_num)
        {
            Point pt = Array.Find(Data.ToArray(), ele => Math.Abs(ele.mass - mass_num) < Constants.tolerance);
            if (pt != null)
            {
                return pt.intensity;
            }
            return -1;
        }

        public void AddData(double mass, double intensity)
        {
            Point pt = new Point();
            pt.mass = mass;
            pt.intensity = intensity;
            Data.Add(pt);
        }

        public double GetMass(int index)
        {
            return Data[index].mass;
        }

        public double GetIntensity(int index)
        {
            return Data[index].intensity;
        }

        public int calculateCharge(double mz)
        {
            for (int i = Charge; i >= 1; i--) // start from highest charge and go to lowest, because say the thing was actually charge 2, and we tested
            {                                 // for charge 1 first. It would register as charge 1 bc it would find a peak every 1 or so.
                bool isotope1 = Array.Exists(Data.ToArray(), element => element.mass >= mz + Constants.isotopeSpacing / i - Constants.tolerance && element.mass <= mz + Constants.isotopeSpacing / i + Constants.tolerance);
                //bool isotope2 = Array.Exists(Data.ToArray(), element => element.mass >= mz + 2*Constants.isotopeSpacing / i - Constants.tolerance && element.mass <= mz + 2*Constants.isotopeSpacing / i + Constants.tolerance);
                if (isotope1)// && isotope2)
                {
                    return i;
                }
            }
            return -1;
        }

        public bool hasGlycanOxoniumIon()
        {
            bool inMassRange1 = Array.Exists(Data.ToArray(), element => element.mass >= Constants.glycan_masses["HexNAc"] + Constants.adduct - Constants.tolerance && element.mass <= Constants.glycan_masses["HexNAc"] + Constants.adduct + Constants.tolerance);
            bool inMassRange2 = Array.Exists(Data.ToArray(), element => element.mass >= Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct - Constants.tolerance && element.mass <= Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct + Constants.tolerance);

            return inMassRange1 || inMassRange2;
        }

        public void sortPeaksByIntensity()
        {
            sortedData = new List<Point>(Data); // contains only most intense peaks

            for (int i = 0; i < sortedData.Count; i++)
            {
                if (sortedData[i].mass <= 600)
                {
                    sortedData.RemoveAt(i);
                    i--;
                }
            }

            sortedData.Sort(new ComparerByIntensity());

            if (sortedData.Count > 25)
            {
                for (int i = 24; i < sortedData.Count;)// for (int i = sortedData.Count; i >=50; i--) may/may not be faster
                {
                    sortedData.RemoveAt(i);
                }
            }
        }

        public Point findGlycopep()
        {
            for (int i = 0; i < sortedData.Count; i++)
            {
                plusHexnac = sortedData[i];  // hypothesized pep+HexNAc peak
                z = calculateCharge(plusHexnac.mass);
                plusHexnac = Utilities.monoisotope(plusHexnac.mass, z, Data);

                mz = Utilities.monoisotope(plusHexnac.mass - Constants.glycan_masses["HexNAc"]/z, z, Data);
                plus2Hexnac = Utilities.monoisotope(mz.mass + (2 * Constants.glycan_masses["HexNAc"])/z, z, Data);
                plusHex2Hexnac = Utilities.monoisotope(mz.mass + (2 * Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"])/z, z, Data);

                bool mz_exists = Array.Exists(Data.ToArray(), element => element.mass >= mz.mass - Constants.tolerance && element.mass <= mz.mass + Constants.tolerance);
                bool mz_2hexNac_exists = Array.Exists(Data.ToArray(), element => element.mass >= plus2Hexnac.mass - Constants.tolerance && element.mass <= plus2Hexnac.mass + Constants.tolerance);
                bool mz_2hexNac_hex_exists = Array.Exists(Data.ToArray(), element => element.mass >= plusHex2Hexnac.mass - Constants.tolerance && element.mass <= plusHex2Hexnac.mass + Constants.tolerance);

                if (mz_exists && mz_2hexNac_exists && mz_2hexNac_hex_exists)
                {
                    return mz;
                }
            }

            Point dun = new Point();
            dun.mass = -1;
            return dun;
        }

        public Point returnHexnac()
        {
            return plusHexnac;
        }

        public Point return2Hexnac()
        {
            return plus2Hexnac;// - _diff;
        }

        public Point returnHex2Hexnac()
        {
            return plusHex2Hexnac;// - _diff;
        }

        public int returnBarePepCharge()
        {
            return z;
        }
    }
}