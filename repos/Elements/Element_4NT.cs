using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elements
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
            this.NumDOFPNode = 2;
            this.NumNodes = 4;
            this.NodalLocations = new double[8];
            this.NodalDisplacements = new double[8];
            this.InternalForce = new double[8];
            this.DMatrix = new(3, 3);
            this.KMatrix = new(8, 8);
        }
    }
}