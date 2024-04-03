namespace FEMAssembly
{
    /// <summary>
    /// 4-Noded Triangular Element
    /// </summary>
    public class Element_4NT : Elements
    {
        // Constructor
        public Element_4NT()
        {
            this.Type = "4NT";
            this.NumDim = 2;
            this.NDOFPNode = 2;
            this.NumNodes = 4;
            this.TotalDOF = this.NDOFPNode * this.NumNodes;
            this.NumIPs = 3;
            this.NodalLocations = new double[this.NDOFPNode * this.NumNodes];
            this.NodalDisplacements = new double[this.NDOFPNode * this.NumNodes];
            this.InternalForce = new double[this.NDOFPNode * this.NumNodes];
            this.ForceVector = new double[this.NDOFPNode * this.NumNodes];
            this.KMatrix = new double[this.NDOFPNode * this.NumNodes, this.NDOFPNode * this.NumNodes];
        }
    }
}