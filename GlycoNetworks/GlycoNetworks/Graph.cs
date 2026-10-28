// © 2026 Alyssa Tang <alyssayztang@gmail.com>

namespace GlycoNetworks
{
    internal class Graph<T>
    {
        public int numNodes;
        public int numEdges;
        public List<T> nodes;
        public List<Edge<T>> edges;

        // merges this graph with another graph, as node1 from this graph is connected to node2 from the otherGraph.
        public virtual void mergeGraph(Graph<T> otherGraph, T node1, T node2, string labelname)
        {
            this.numNodes += otherGraph.numNodes;
            this.numEdges += otherGraph.numEdges;
            // merged node list from otherGraph
            for (int i = 0; i < otherGraph.numNodes; i++)
            {
                this.nodes.Add(otherGraph.nodes[i]);
            }
            // merged edge list from otherGraph 
            for (int i = 0; i < otherGraph.numEdges; i++)
            {
                this.edges.Add(otherGraph.edges[i]);
            }
            // adding edges connecting node1 to node2 + vice versa
            this.addEdge(node1, node2, labelname);
        }

        public virtual void addEdge(T node1, T node2, string labelname)
        {
            Edge<T> edge = new Edge<T>();
            edge.source = node1;
            edge.target = node2;
            edge.glycDiff = labelname;
            this.edges.Add(edge);
        }
    }

    internal class GlycoGraph : Graph<Glycopep>
    {
        public override void addEdge(Glycopep node1, Glycopep node2, string labelname)
        {
            Glycopep biggerGlycan, smallerGlycan;
            if (node1.glycanMass > node2.glycanMass)
            {
                biggerGlycan = node1;
                smallerGlycan = node2;
            }
            else
            {
                biggerGlycan = node2;
                smallerGlycan = node1;
            }

            Edge<Glycopep> edge1 = new Edge<Glycopep>();
            edge1.source = biggerGlycan;
            edge1.target = smallerGlycan;
            edge1.glycDiff = labelname;
            edge1.increasing = -1;

            Edge<Glycopep> edge2 = new Edge<Glycopep>();
            edge2.source = smallerGlycan;
            edge2.target = biggerGlycan;
            edge2.glycDiff = labelname;
            edge2.increasing = 1;

            this.edges.Add(edge1);
            this.edges.Add(edge2);
        }
    }

    // previously used in the deprecated method without diffValue buckets 
    internal class TupleGraph : Graph<Tuple<double, int>>
    {
        public Tuple<double, int> largestNode;

        public void mergeGraph(TupleGraph otherGraph, Tuple<double, int> node1, Tuple<double, int> node2, string labelname)
        {
            base.mergeGraph(otherGraph, node1, node2, labelname);

            if (this.largestNode.Item2 < otherGraph.largestNode.Item2)
            {
                this.largestNode = otherGraph.largestNode;
            }
            else if (this.largestNode.Item2 == otherGraph.largestNode.Item2 && Math.Abs(this.largestNode.Item1) > Math.Abs(otherGraph.largestNode.Item1))
            {
                this.largestNode = otherGraph.largestNode;
            }
        }
    }

    internal class DiffValueGraph : Graph<Tuple<double, Bucket<DiffValue>>>
    {

        public Tuple<double, Bucket<DiffValue>> largestNode;

        public void mergeGraph(DiffValueGraph otherGraph, Tuple<double, Bucket<DiffValue>> node1, Tuple<double, Bucket<DiffValue>> node2, string labelname)
        {
            base.mergeGraph(otherGraph, node1, node2, labelname);

            if (this.largestNode.Item2.items.Count < otherGraph.largestNode.Item2.items.Count)
            {
                this.largestNode = otherGraph.largestNode;
            }
            else if (this.largestNode.Item2.items.Count == otherGraph.largestNode.Item2.items.Count && Math.Abs(this.largestNode.Item2.items.Count) > Math.Abs(otherGraph.largestNode.Item2.items.Count))
            {
                this.largestNode = otherGraph.largestNode;
            }
        }
    }

    internal class MergedBlockGraph : Graph<MergedBlock>
    {
        public MergedBlock largestNode;

        public void mergeGraph(MergedBlockGraph otherGraph, MergedBlock node1, MergedBlock node2, string labelname)
        {
            base.mergeGraph(otherGraph, node1, node2, labelname);

            if (this.largestNode.confidence < otherGraph.largestNode.confidence)
            {
                this.largestNode = otherGraph.largestNode;
            }
        }
    }
}