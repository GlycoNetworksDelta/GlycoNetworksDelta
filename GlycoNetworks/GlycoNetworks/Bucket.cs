// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;

namespace GlycoNetworks
{
    internal class Bucket<T>
    {
        public HashSet<T> items;
        public List<T> itemsList;
        public T medianItem;

        public Bucket()
        {
            items = new HashSet<T>();
            itemsList = new List<T>();
        }

        public Bucket(T item)
        {
            items = new HashSet<T>() { item };
            itemsList = new List<T>() { item };
        }

        public int calculateMedianIndex()
        {
            int len = itemsList.Count;
            int medianIndex = len / 2;
            return medianIndex;
        }

        public virtual void calculateMedianItem()
        {
            int medianIndex = calculateMedianIndex();
            itemsList = itemsList.Order().ToList();

            medianItem = itemsList[medianIndex];
        }
    }

    internal class BlockBucket : Bucket<Block>
    {
        public List<double> mergedRealMassValues;
        public List<realMassDifferences> mergedRealMassDifferences;
        public double mergedConfidence;
        public HashSet<int> oxoniumIonMasses; 
        public DescriptiveStatistics statistics;
        public double marginOfError;
        public double marginOfError2;
        public double medianValue;

        //debugging
        public double numSpectraTotal;

        public BlockBucket()
        {
            items = new HashSet<Block>();
            itemsList = new List<Block>();
            mergedRealMassValues = new List<double>();
            mergedRealMassDifferences = new List<realMassDifferences>();
            mergedConfidence = 0.0;
            oxoniumIonMasses = new HashSet<int>();
        }

        public BlockBucket(Block item)
        {
            items = new HashSet<Block>() { item };
            itemsList = new List<Block>() { item };
            mergedRealMassValues = new List<double>();
            mergedRealMassDifferences = new List<realMassDifferences>();
            mergedConfidence = 0.0;
            oxoniumIonMasses = new HashSet<int>();
        }

        public List<double> massValues()
        {
            var massValues = new List<double>();
            foreach (Block block in itemsList)
            {
                massValues.Add(block.mass);
            }
            return massValues;
        }
    }
}
