// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>
// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal static class Constants
    {
        public static double tolerance = 0.03;

        public static Dictionary<String, double> glycan_masses = new Dictionary<String, double>
        {
            ["Hex"] = 162.052823,
            ["HexNAc"] = 203.079373,
            ["Fuc"] = 146.057909,
            ["NeuAc"] = 291.095417,
            ["NeuGc"] = 307.090331,
            ["Pent"] = 150.05282342 - 18.010565,
            ["Phospho"] = 79.966331,
            ["Na"] = 21.981943,
            ["Acetyl"] = 42.010565,
            ["Sulfo"] = 79.956815,
            ["pC"] = 165.055480,
            ["Hex-HexNAc"] = 162.052823 + 203.079373  // Special case of the most common disaccharide delta
        };

        public static Dictionary<String, double> glycanSpecificIons = new Dictionary<String, double>
        {
            ["Hex"] = 162.052823 + adduct,
            ["HexNAc"] = 203.079373 + adduct,
            ["Hex-HexNAc"] = 162.052823 + 203.079373 + adduct,
            ["NeuAc"] = 291.095417 + adduct,
            ["NeuGc"] = 307.090331 + adduct
        };

        public const double isotopeSpacing = 1.003;
        public const double adduct = 1.007825035 - 0.00054858;
        public const double water = 18.010565;

        public static double weight = 1; //the smaller the weight, the closer confidence is to deltaHi
        // 0 is when its deltaHits 
        // anything below 0 it becomes like.... more weight for the peaks now. 
        public static int verbosityLevel = 1;
    }
}
