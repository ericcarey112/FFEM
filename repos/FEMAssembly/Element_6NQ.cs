namespace FEMAssembly
{
    /// <summary>
    /// 6-Noded Quadrilateral Element
    /// </summary>
    public class Element_6NQ : Elements
    {
        // Constructor
        public Element_6NQ()
        {
            this.Type = "6NQ";
            this.NumDim = 2;
            this.NDOFPNode = 2;
            this.NumNodes = 6;
            this.TotalDOF = this.NDOFPNode * this.NumNodes;
            this.NumIPs = 9;
            this.NodalLocations = new double[this.NDOFPNode * this.NumNodes];
            this.NodalDisplacements = new double[this.NDOFPNode * this.NumNodes];
            this.InternalForce = new double[this.NDOFPNode * this.NumNodes];
            this.ForceVector = new double[this.NDOFPNode * this.NumNodes];
            this.KMatrix = new double[this.NDOFPNode * this.NumNodes, this.NDOFPNode * this.NumNodes];
        }
    }
}