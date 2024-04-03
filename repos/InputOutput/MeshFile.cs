//using FEMAssembly;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEMAssembly
{
    /// <summary>
    /// Mesh text file generated in MATLAB
    /// </summary>
    public class MeshFile
    {
        // Properties
        public string FileName { get; set; }
        public string FileDirectory { get; set; }

        // Constructor
        public MeshFile(string filename, string directory)
        {
            this.FileName = filename;
            this.FileDirectory = directory;
        }

        // Methods
        /// <summary>
        /// Read in the nodal locations, connectivity, element sets, periodic node pairs, and pinned nodes from text file
        /// </summary>
        public void ReadMeshFile(Assembly assembly)
        {
            string FilePath = Path.Combine(this.FileDirectory, this.FileName);

            // Make sure file exists
            if (File.Exists(FilePath))
            {
                using StreamReader reader = new(FilePath);
                string[] lines = File.ReadAllLines(FilePath);
                bool NodeFlag = false;
                bool ConnectivityFlag = false;
                bool NumElementsFlag = false;
                bool Fiber4NTFlag = false;
                bool Matrix6NTFlag = false;
                bool Matrix8NQFlag = false;
                bool RightEdgeFlag = false;
                bool TopEdgeFlag = false;
                bool NodePairFlag = false;
                bool PinnedNodesFlag = false;
                bool RVELengthFlag = false;
                foreach (string line in lines)
                {
                    // Determine current file sections
                    if (line == "**NODES") { NodeFlag = true; continue; }
                    else if (line == "**CONNECTIVITY") { NodeFlag = false; ConnectivityFlag = true; continue; }
                    else if (line == "**NUMBER OF ELEMENTS") { ConnectivityFlag = false; NumElementsFlag = true; continue; }
                    else if (line == "**FIBER ELEMENTS 4NT") { NumElementsFlag = false; Fiber4NTFlag = true; continue; }
                    else if (line == "**MATRIX ELEMENTS 6NT") { Fiber4NTFlag = false; Matrix6NTFlag = true; continue; }
                    else if (line == "**MATRIX ELEMENTS 8NQ") { Matrix6NTFlag = false; Matrix8NQFlag = true; continue; }
                    else if (line == "**RIGHT EDGE") { Matrix8NQFlag = false; RightEdgeFlag = true; continue; }
                    else if (line == "**TOP EDGE") { RightEdgeFlag = false; TopEdgeFlag = true; continue; }
                    else if (line == "**PERIODIC NODE PAIRS") { TopEdgeFlag = false; NodePairFlag = true; continue; }
                    else if (line == "**PINNED NODES") { NodePairFlag = false; PinnedNodesFlag = true; continue; }
                    else if (line == "**RVE LENGTH") { PinnedNodesFlag = false; RVELengthFlag = true; continue; }
                    else if (line == "**END")
                    {
                        assembly.NumberOfNodes = assembly.NodalLocations.Count;
                        assembly.TotalDOF = assembly.NDOFPNode * assembly.NumberOfNodes;
                        assembly.GlobalQ = new double[assembly.TotalDOF];
                        return;
                    }

                    // Read in mesh data
                    if (NodeFlag) { ReadGlobalNodes(line, assembly); }
                    if (ConnectivityFlag) { ReadConnecivity(line, assembly); }
                    if (NumElementsFlag) { ReadNumberOfElements(line, assembly); }
                    if (Fiber4NTFlag) { ReadElementType(line, assembly, "4NT"); }
                    if (Matrix6NTFlag) { ReadElementType(line, assembly, "6NT"); }
                    if (Matrix8NQFlag) { ReadElementType(line, assembly, "8NQ"); }
                    if (RightEdgeFlag) { ReadEdgeNodes(line, assembly, "Right"); }
                    if (TopEdgeFlag) { ReadEdgeNodes(line, assembly, "Top"); }
                    if (NodePairFlag) { ReadNodePairArray(line, assembly); }
                    if (PinnedNodesFlag) { ReadPinnedNodes(line, assembly); }
                    if (RVELengthFlag) { ReadRVELength(line, assembly); }
                }
            }
            else
            {
                throw new Exception("Could not find mesh file. Check file name and directory");
            }
        }

        /// <summary>
        /// Read in global node locations
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadGlobalNodes(string line, Assembly assembly)
        {
            string[] nodeString = line.Split(' ');
            double[] nodeDouble = new double[3];
            for (int i = 0; i < nodeDouble.Length; i++)
            {
                if (double.TryParse(nodeString[i], out double coordinate))
                {
                    nodeDouble[i] = coordinate;
                }
                else
                {
                    throw new Exception("Error in reading global node coordinates in InputOutput --> ReadGlobalNodes");
                }
            }
            assembly.NodalLocations.Add(nodeDouble);
        }

        /// <summary>
        /// Read in element connectivity
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadConnecivity(string line, Assembly assembly)
        {
            string[] nodeString = line.Split(' ');
            int[] nodeInt = new int[nodeString.Length];
            for (int i = 0; i < nodeString.Length; i++)
            {
                if (int.TryParse(nodeString[i], out int node))
                {
                    nodeInt[i] = node;
                }
                else
                {
                    throw new Exception("Error in reading connectivity in InputOutput --> ReadConnectivity");
                }
            }
            assembly.Connectivity.Add(nodeInt);
        }

        /// <summary>
        /// Read in number of elements
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadNumberOfElements(string line, Assembly assembly)
        {
            if (int.TryParse(line, out int NumElements))
            {
                assembly.NumberOfElements = NumElements;
                assembly.ElementList = new List<Elements>(NumElements);
                // Add place holders in list so we can directly assign elements based on their number
                Elements element = new Element_4NQ();
                for (int i = 0; i < NumElements; i++) { assembly.ElementList.Add(element); }
            }
            else
            {
                throw new Exception("Error in reading number of elements in InputOutput --> ReadNumberOfElementsNodes");
            }
        }

        /// <summary>
        /// Create an element object based on type and store in assembly.ElementList
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <param name="type"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadElementType(string line, Assembly assembly, string type)
        {
            int ElemNum;
            Elements element;

            // 4-Noded Triangles (Fibers)
            if (type == "4NT")
            {
                if (int.TryParse(line, out ElemNum))
                {
                    element = new Element_4NT();
                }
                else
                {
                    throw new Exception("Error in assigning 4NT element type in InputOutput --> ReadElementTypes()");
                }
            }

            // 6-Noded Triangles (Matrix)
            else if (type == "6NT")
            {
                if (int.TryParse(line, out ElemNum))
                {
                    element = new Element_6NT();
                }
                else
                {
                    throw new Exception("Error in assigning 6NT element type in InputOutput --> ReadElementTypes()");
                }
            }

            // 8-Noded Quadrilaterals (Matrix between fibers)
            else if (type == "8NQ")
            {
                if (int.TryParse(line, out ElemNum))
                {
                    element = new Element_8NQ();
                }
                else
                {
                    throw new Exception("Error in assigning 8NQ element type in InputOutput --> ReadElementTypes()");
                }
            }
            else { throw new Exception("Error in element type name in InputOutput --> ReadElementTypes()"); }
            AssignNodalLocationsToElement(assembly, element, ElemNum);
            assembly.ElementList[ElemNum - 1] = element;
            assembly.ElementList[ElemNum - 1].ElementNumber = ElemNum;
        }

        /// <summary>
        /// Assign element nodal locations based on connectivity and global node locations
        /// </summary>
        private static void AssignNodalLocationsToElement(Assembly assembly, Elements element, int ElemNum)
        {
            int count = 0;
            for (int i = 0; i < element.NumNodes; i++)
            {
                int GlobalNodeNum = assembly.Connectivity[ElemNum - 1][i];
                for (int j = 0; j < element.NumDim; j++)
                {
                    element.NodalLocations[count] = assembly.NodalLocations[GlobalNodeNum - 1][j];
                    count++;
                }
            }
        }

        /// <summary>
        /// Assign edge nodes (Top and Right edge nodes) for use in homogenization
        /// </summary>
        private static void ReadEdgeNodes(string line, Assembly assembly, string Edge)
        {
            // Nodes along right edge of RVE
            if (Edge == "Right")
            {
                if (int.TryParse(line, out int NodeNum))
                {
                    assembly.RightEdgeNodes.Add(NodeNum);
                }
                else
                {
                    throw new Exception("Error in assigning right edge nodes in InputOutput --> ReadElementTypes()");
                }
            }

            // Nodes along top edge
            else if (Edge == "Top")
            {
                if (int.TryParse(line, out int NodeNum))
                {
                    assembly.TopEdgeNodes.Add(NodeNum);
                }
                else
                {
                    throw new Exception("Error in assigning top edge nodes in InputOutput --> ReadElementTypes()");
                }
            }
        }

        /// <summary>
        /// Read periodic node pair array
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadNodePairArray(string line, Assembly assembly)
        {
            string[] pairString = line.Split(' ');
            int[] nodeInt = new int[pairString.Length];
            for (int i = 0; i < pairString.Length; i++)
            {
                if (int.TryParse(pairString[i], out int node))
                {
                    nodeInt[i] = node;
                }
                else
                {
                    throw new Exception("Error in reading node pair array in InputOutput --> ReadNodePairArray");
                }
            }
            assembly.NodePairArray.Add(nodeInt);
        }

        /// <summary>
        /// Read in pinned nodes for BCs
        /// </summary>
        private static void ReadPinnedNodes(string line, Assembly assembly)
        {
            if (int.TryParse(line, out int node))
            {
                assembly.PinnedNode = node;
            }
            else
            {
                throw new Exception("Error in reading node pair array in InputOutput --> ReadNodePairArray");
            }
        }

        /// <summary>
        /// Read in RVE Lengths (2 and 3 direction)
        /// </summary>
        /// <param name="line"></param>
        /// <param name="assembly"></param>
        /// <exception cref="Exception"></exception>
        private static void ReadRVELength(string line, Assembly assembly)
        {
            string[] lengthString = line.Split(' ');
            for (int i = 0; i < lengthString.Length; i++)
            {
                if (double.TryParse(lengthString[i], out double length))
                {
                    assembly.RVELength[i] = length;
                }
                else
                {
                    throw new Exception("Error in reading RVE length in InputOutput --> ReadNodePairArray");
                }
            }
        }
    }
}
