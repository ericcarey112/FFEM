namespace FEMAssembly
{
    /// <summary>
    /// 3-Noded Triangular Element
    /// </summary>
    public class Element_3NT : Elements
    {
        // Constructor
        public Element_3NT()
        {
            this.Type = "3NT";
            this.NumDim = 2;
            this.NDOFPNode = 2;
            this.NumNodes = 3;
            this.TotalDOF = this.NDOFPNode * this.NumNodes;
            this.NumIPs = 1;
            this.NodalLocations = new double[this.NDOFPNode * this.NumNodes];
            this.NodalDisplacements = new double[this.NDOFPNode * this.NumNodes];
            this.InternalForce = new double[this.NDOFPNode * this.NumNodes];
            this.ForceVector = new double[this.NDOFPNode * this.NumNodes];
            this.KMatrix = new double[this.NDOFPNode * this.NumNodes, this.NDOFPNode * this.NumNodes];
        }
    }
}