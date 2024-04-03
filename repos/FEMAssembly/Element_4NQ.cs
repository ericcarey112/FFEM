namespace FEMAssembly
{
    /// <summary>
    /// 4-Noded Quadrilateral Element
    /// </summary>
    public class Element_4NQ : Elements
    {
        // Constructor
        public Element_4NQ()
        {
            this.Type = "4NQ";
            this.NumDim = 2;
            this.NDOFPNode = 2;
            this.NumNodes = 4;
            this.TotalDOF = this.NDOFPNode * this.NumNodes;
            this.NumIPs = 4;
            this.NodalLocations = new double[this.NDOFPNode * this.NumNodes];
            this.NodalDisplacements = new double[this.NDOFPNode * this.NumNodes];
            this.InternalForce = new double[this.NDOFPNode * this.NumNodes];
            this.ForceVector = new double[this.NDOFPNode * this.NumNodes];
            this.KMatrix = new double[this.NDOFPNode * this.NumNodes, this.NDOFPNode * this.NumNodes];
        }
    }
}
