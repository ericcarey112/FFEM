﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elements
{
    /// <summary>
    /// 6-Noded Triangular Element
    /// </summary>
    public class Element_6NT : Elements
    {
        // Constructor
        public Element_6NT()
        {
            this.Type = "6NT";
            this.NumDim = 2;
            this.NumDOFPNode = 2;
            this.NumNodes = 6;
            this.NodalLocations = new double[12];
            this.NodalDisplacements = new double[12];
            this.InternalForce = new double[12];
            this.DMatrix = new(3, 3);
            this.KMatrix = new(12, 12);
        }
    }
}