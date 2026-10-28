// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>
// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal class GlycanComposition : Dictionary<String, int>
    {
        public String toString()
        {
            String strComp = "";
            for (int i = 0; i < Keys.Count; i++)
            {
                if (this[Keys.ToList()[i]] != 0)
                {
                    string key = Keys.ToList()[i];
                    int count = this[Keys.ToList()[i]];
                    double keyAsMass;
                    if (Double.TryParse(key, out keyAsMass))
                    {
                        strComp += ' ' + (count * keyAsMass).ToString();
                    }
                    else
                    {
                        strComp += key + "(" + count + ")";
                    }
                }
            }
            return strComp;
        }

        public bool isValid()
        {
            foreach (KeyValuePair<String, int> kvp in this)
            {
                if (kvp.Value < 0) return false;
            }
            return true;
        }

        public static GlycanComposition toComp(String input)
        {
            String strComp = input;
            GlycanComposition comp = new GlycanComposition();
            char[] delimiters = { '(', ')' };
            String[] listComp = strComp.Split(delimiters);
            String key = "";
            int val = 0;
            for (int i = 0; i < listComp.Length; i++)
            {
                if (i % 2 == 0)
                {
                    key = listComp[i];
                }
                else
                {
                    val = Convert.ToInt32(listComp[i]);
                    comp.Add(key, val);
                }

                //int val = Convert.ToInt32(strComp.Split("(")[1].Split(")")[0]);
                //Console.Write(key);
                //Console.WriteLine(val);
                //comp.Add(key, val);
                //strComp = strComp.Split(key + "(" + val.ToString() + ")")[1];
                //Console.WriteLine(strComp);
                //Console.Write(key);
                //Console.WriteLine(val);
                //comp.Add(key, val);
            }
            return comp;
        }

        public static bool Contains(List<GlycanComposition> list, GlycanComposition x)
        {
            foreach (GlycanComposition listItem in list)
            {
                if (x.toString() == listItem.toString()) return true;
            }
            return false;
        }

        public static bool ContainsMassMatch(List<GlycanComposition> list, GlycanComposition x)
        {
            var compareGlycan = new CompareGlycan();
            foreach (GlycanComposition listItem in list)
            {
                if (Math.Abs(compareGlycan.calcMass(x) - compareGlycan.calcMass(listItem)) < 0.00001) return true;
            }
            return false;
        }
    }
}
