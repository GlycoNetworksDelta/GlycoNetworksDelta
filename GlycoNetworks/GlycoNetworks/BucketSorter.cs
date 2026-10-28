// © 2026 Alyssa Tang <alyssayztang@gmail.com>

using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;

namespace GlycoNetworks
{
    internal class BucketSorter<T> 
    {
        public Dictionary<double, Bucket<T>> buckets;
        public int numBuckets;
    }

    internal class DiffValueBuckets : BucketSorter<DiffValue>
    {
        public DiffValueBuckets(List<DiffValue> items)
        {
            buckets = new Dictionary<double, Bucket<DiffValue>>();
            this.addingItems(items);
            numBuckets = buckets.Count;
        }

        public DiffValueBuckets(List<DiffValue> items, List<double> starterKeys)
        {
            buckets = new Dictionary<double, Bucket<DiffValue>>();
            foreach (double starterKey in starterKeys)
            {
                Boolean keyAdded = false;
                foreach (double key in buckets.Keys)
                {
                    if (starterKey - Constants.tolerance < key && key < starterKey + Constants.tolerance)
                    {
                        keyAdded = true;
                    }
                }
                if (keyAdded == false) {
                    buckets[starterKey] = new Bucket<DiffValue>();
                }
            }
            this.addingItems(items);
            numBuckets = buckets.Count;
        }

        public void addingItems(List<DiffValue> items)
        {
            foreach (DiffValue item in items)
            {
                Boolean itemAdded = false;
                foreach (double key in buckets.Keys)
                {
                    if (item.diff - Constants.tolerance < key && key < item.diff + Constants.tolerance)
                    {
                        buckets[key].items.Add(item);
                        buckets[key].itemsList.Add(item);
                        itemAdded = true;
                    }
                }
                if (itemAdded == false)
                {
                    buckets[item.diff] = new Bucket<DiffValue>(item);
                }
            }
        }

        public List<double> returnMedians()
        {
            var medians = new List<double>();

            foreach (double key in buckets.Keys)
            {
                List<double> diffs = new List<double>();
                foreach (DiffValue item in  buckets[key].itemsList)
                {
                    diffs.Add(item.diff);
                }
                double medianMass = MathNet.Numerics.Statistics.Statistics.Median(diffs);
                medians.Add(medianMass);
            }

            return medians;
        }
    }

    internal class BlockBuckets : BucketSorter<Block>
    {

        public Dictionary<double, BlockBucket> blockBuckets;

        public BlockBuckets(List<Block> items)
        {
            blockBuckets = new Dictionary<double, BlockBucket>();
            this.addingItems(items);
            numBuckets = blockBuckets.Count;
        }

        public BlockBuckets(List<Block> items, List<double> starterKeys)
        {
            blockBuckets = new Dictionary<double, BlockBucket>();
            foreach (double starterKey in starterKeys)
            {
                Boolean keyAdded = false;
                foreach (double key in blockBuckets.Keys)
                {
                    if (starterKey - Constants.tolerance < key && key < starterKey + Constants.tolerance)
                    {
                        keyAdded = true;
                    }
                }
                if (keyAdded == false) {
                    blockBuckets[starterKey] = new BlockBucket();
                }
            }
            this.addingItems(items);
            numBuckets = blockBuckets.Count;
        }

        public void addingItems(List<Block> items)
        {
            foreach (Block item in items)
            {
                Boolean itemAdded = false;
                foreach (double key in blockBuckets.Keys)
                {
                    if (item.mass - Constants.tolerance < key && key < item.mass + Constants.tolerance)
                    {
                        blockBuckets[key].items.Add(item);
                        blockBuckets[key].itemsList.Add(item);
                        itemAdded = true;
                    }
                }
                if (itemAdded == false)
                {
                    blockBuckets[item.mass] = new BlockBucket(item);
                }
            }
        }

        public List<double> returnMedians()
        {
            var medians = new List<double>();

            foreach (double key in blockBuckets.Keys)
            {
                double medianMass = MathNet.Numerics.Statistics.Statistics.Median(blockBuckets[key].massValues());
                medians.Add(medianMass);
            }

            return medians;
        }
    }
}
