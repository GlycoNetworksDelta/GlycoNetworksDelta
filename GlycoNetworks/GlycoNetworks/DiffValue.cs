// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal class DiffValue
    {
        public double diff; // diff is graphJ node - graphI node 

        public Glycopep nodeGraphI;
        public Glycopep nodeGraphJ;
        public GlycoGraph graphI;
        public GlycoGraph graphJ;

        public DiffValue(double diff, GlycoGraph graphI, GlycoGraph graphJ, Glycopep nodeGraphI, Glycopep nodeGraphJ)
        {
            this.diff = diff;
            this.nodeGraphI = nodeGraphI;
            this.nodeGraphJ = nodeGraphJ;
            this.graphI = graphI;
            this.graphJ = graphJ;
        }
    }
}
