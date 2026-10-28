// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.Distributions;

namespace GlycoNetworks
{
    internal class Block
    {
        public double mass;

        public int diffValueHits; // number of hits from the diff value graphs
        public DiffValueGraph diffValueGraph;
        // public GlycoGraph smallerGlycograph; //aka the control graph
        // public GlycoGraph largerGlycograph; // graph2 is the one w the builidng block 
        //better naming scheme needed ;-; 
        public HashSet<Glycopep> smallerGlycographRelevantNodes;
        public HashSet<Glycopep> largerGlycographRelevantNodes;

        public int deltaHits; // number of hits from examining delta peaks in spectra
        public List<Tuple<MS2, Point>> controlDeltaPeaks; //pair of spectra (MS2) and point from spectra, for delta peaks
        public List<Tuple<MS2, Point>> controlDeltaBlockPeaks;
        public List<Tuple<MS2, Point>> deltaPeaks;
        public List<Tuple<MS2, Point>> deltaBlockPeaks; // like deltaPeaks, but for delta + this building block peaks 

        public double controlDeltaConfidence;
        public double controlDeltaBlockConfidence;
        public double deltaConfidence;
        public double deltaBlockConfidence;
        public double confidence; //deltaBlockConfidence - controlDeltaBlockConfidence 

        // what do i want for my metadata?
        // want: the specific value differences
        // the spectra 
        public List<realMassDifferences> realMassDifferences; 
        public List<double> realMassValues;
        public double medianRealMassValue;
        //public DescriptiveStatistics statistics;
        //public double marginOfError;

        // debugging
        public double numSpectra;
        public List<MS2> spectra;
        public List<int> indices;
        public Dictionary<MS2, List<double>> RTToMass;

    }

    internal class MergedBlock : Block
    {
        public Tuple<double, BlockBucket> keyValuePair; // tho technically, the current mass is the key in the bucketsorter 
        public BlockBucket blockBucket; 
    }

    internal class  realMassDifferences
    {
        public double diff;
        public double smallerValue;
        public double largerValue;
        public MS2 spectra; 

        public realMassDifferences(double diff, double smallerValue, double largeerValue, MS2 spectra)
        {
            this.diff = diff;
            this.smallerValue = smallerValue;
            this.largerValue = largeerValue;
            this.spectra = spectra;
        }
    }
}
