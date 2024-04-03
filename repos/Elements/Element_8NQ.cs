using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elements
{
    /// <summary>
    /// 8-Noded Quadrilateral Element
    /// </summary>
    public class Element_8NQ : Elements
    {
        // Constructor
        public Element_8NQ()
        {
            this.Type = "8NQ";
            this.NumDim = 2;
            this.NumDOFPNode = 2;
            this.NumNodes = 8;
            this.NodalLocations = new double[16];
            this.NodalDisplacements = new double[16];
            this.InternalForce = new double[16];
            this.DMatrix = new(3, 3);
            this.KMatrix = new(32, 32);
        }
    }
}