using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elements
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
            this.NumDOFPNode = 2;
            this.NumNodes = 3;
            this.NodalLocations = new double[6];
            this.NodalDisplacements = new double[6];
            this.InternalForce = new double[6];
            this.DMatrix = new(3, 3);
            this.KMatrix = new(6, 6);
        }
    }
}