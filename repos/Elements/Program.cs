using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

namespace Elements
{
    class Program
    {
        public static void Main(string[] args)
        {
            // Construct Element
            Elements Elem1 = new Element_3NT
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 1
            };

            Elements Elem2 = new Element_4NT
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 3
            };

            Elements Elem3 = new Element_6NT
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 7
            };

            Elements Elem4 = new Element_4NQ
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 4
            };

            Elements Elem5 = new Element_6NQ
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 4
            };

            Elements Elem6 = new Element_8NQ
            {
                // Material Constants and IPs
                E2 = 70000,
                nu23 = 0.3,
                NumIPs = 9
            };

            // Nodal Locations and Displacements
            double[] Locations1 = [0,0, 1,0, 0,1];
            double[] Locations2 = [0,1, 0,0, 0.5,0, 1,0];
            double[] Locations3 = [0,0, 0.5,0, 1,0, 0.5,0.5, 0,1, 0,0.5];
            double[] Locations4 = [0,0, 1,0, 1,1, 0,1];
            double[] Locations5 = [0,0, 0.5,0, 1,0, 1,1, 0.5,1, 0,1];
            double[] Locations6 = [0,0, 0.5,0, 1,0, 1,0.5, 1,1, 0.5,1, 0,1, 0,0.5];
            double[] Displacements1 = [0,0, 1,0, 0,0];
            double[] Displacements2 = [0,0, 0,1, 0,1, 0,1];
            double[] Displacements3 = [0,1, 0.5,0, 2,0, 0.5,1, 1,1, 0,0.5];
            double[] Displacements4 = [0, 0, 1, 0, 1, 0, 0, 0];
            double[] Displacements5 = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0];
            double[] Displacements6 = [0,0, 0,0, 1,0, 1,0, 1,0, 0,0, 0,0, 0,0];
            Elem1.NodalLocations = Locations1;
            Elem1.NodalDisplacements = Displacements1;
            Elem2.NodalLocations = Locations2;
            Elem2.NodalDisplacements = Displacements2;
            Elem3.NodalLocations = Locations3;
            Elem3.NodalDisplacements = Displacements3;
            Elem4.NodalLocations = Locations4;
            Elem4.NodalDisplacements = Displacements4;
            Elem5.NodalLocations = Locations5;
            Elem5.NodalDisplacements = Displacements5;
            Elem6.NodalLocations = Locations6;
            Elem6.NodalDisplacements = Displacements6;

            // K Matrix
            Elem1.KMatrix2D(Elem1.NodalDisplacements);
            Elem2.KMatrix2D(Elem2.NodalDisplacements);
            Elem3.KMatrix2D(Elem3.NodalDisplacements);
            Elem4.KMatrix2D(Elem4.NodalDisplacements);
            Elem5.KMatrix2D(Elem5.NodalDisplacements);
            Elem6.KMatrix2D(Elem6.NodalDisplacements);

            Console.WriteLine(Elem6.KMatrix);
            
            Console.ReadKey();
        }
    }

    public class Elements
    {
        // Properties
        public string Type { get; set; }
        public int NumDim { get; set; }
        public int NumNodes { get; set; }
        public int NumDOFPNode { get; set; }
        public int NumIPs {  get; set; }
        public double Thickness { get; set; }
        public string Isotropy { get; set; }
        public double E1 {  get; set; }
        public double E2 { get; set; }
        public double E3 { get; set; }
        public double nu12 {  get; set; }
        public double nu23 { get; set; }
        public double G12 {  get; set; }
        public double G23 { get; set; }
        public double[] NodalLocations { get; set; }
        public double[] NodalDisplacements { get; set; }
        public double[] InternalForce { get; set; }
        public Matrix DMatrix { get; set; }
        public Matrix KMatrix { get; set; }
        public double[] Damage { get; set; }

        // Constructor
        public Elements()
        {
            this.Type = " ";
            this.NumDim = 0;
            this.NumNodes = 0;
            this.NumDOFPNode = 0;
            this.NumIPs = 1;
            this.Thickness = 1.0;
            this.Isotropy = " ";
            this.E1 = 70000.0;
            this.E2 = 70000.0;
            this.E3 = 70000.0;
            this.nu12 = 0.3;
            this.nu23 = 0.3;
            this.G12 = 20000.0;
            this.G23 = 20000.0;
            this.NodalLocations = new double[1];
            this.NodalDisplacements = new double[1];
            this.InternalForce = new double[1];
            this.DMatrix = new Matrix(1, 1);
            this.KMatrix = new Matrix(1, 1);
            this.Damage = new double[2];
        }

        // Methods
        /// <summary>
        /// 2D DMatrix for isotropic plane stress element
        /// </summary>
        public void DMatrixPlaneStressIsotropic()
        {
            double dmg = this.Damage[1]; // Damage variable
            double E = this.E2;          // E  (isotropic)
            double nu = this.nu23;       // nu (isotropic)

            // Calculate D
            double constant = E / (1.0 - nu * nu);
            this.DMatrix[0, 0] = (1.0 - dmg) * constant;
            this.DMatrix[0, 1] = (1.0 - dmg) * nu * constant;
            this.DMatrix[0, 2] = 0.0;
            this.DMatrix[1, 0] = this.DMatrix[0, 1];
            this.DMatrix[1, 1] = this.DMatrix[0, 0];
            this.DMatrix[1, 2] = 0.0;
            this.DMatrix[2, 0] = this.DMatrix[0, 2];
            this.DMatrix[2, 1] = this.DMatrix[1, 2];
            this.DMatrix[2, 2] = (1.0 - dmg) * constant * (1 - nu) / (2.0);
        }

        /// <summary>
        /// 2D DMatrix for isotropic plane strain element
        /// </summary>
        public void DMatrixPlaneStrainIsotropic()
        {
            double dmg = this.Damage[1]; // Damage variable
            double E = this.E2;          // E  (isotropic)
            double nu = this.nu23;       // nu (isotropic)

            // Calculate D
            double constant = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            this.DMatrix[0, 0] = (1.0 - dmg) * constant * (1.0 - nu);
            this.DMatrix[0, 1] = (1.0 - dmg) * constant * nu;
            this.DMatrix[0, 2] = 0.0;
            this.DMatrix[1, 0] = this.DMatrix[0, 1];
            this.DMatrix[1, 1] = this.DMatrix[0, 0];
            this.DMatrix[1, 2] = 0.0;
            this.DMatrix[2, 0] = this.DMatrix[0, 2];
            this.DMatrix[2, 1] = this.DMatrix[1, 2];
            this.DMatrix[2, 2] = (1.0 - dmg) * constant * (1 - 2.0 * nu) / (2.0);
        }

        /// <summary>
        /// Calculate K matrix for 2D element
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        public void KMatrix2D(double[] q)
        {
            // Get weights and IPs
            Matrix W_xi_eta = GetIntegrationPoints(this, NumIPs);

            // Calculate K Matrix
            Matrix TotalKMatrix = new(NumNodes * NumDOFPNode, NumNodes * NumDOFPNode);
            double[] TotalInternalForce = new double[NumNodes * NumDOFPNode];
            for (int i = 0; i < NumIPs; i++)
            {
                double W = W_xi_eta[i, 0]; // Weight
                double xi = W_xi_eta[i, 1]; // xi
                double eta = W_xi_eta[i, 2]; // eta

                // Calculate DMatrix
                this.DMatrixPlaneStrainIsotropic();

                // Calculate Jacobian
                Matrix Jacobian = Elements.Jacobian2D(this, xi, eta);
                double detJ = Jacobian.Det();

                // Calculate B Matrix:
                Matrix BMatrix = Elements.BMatrix2D(this, xi, eta);
                Matrix BMatrixTranspose = Matrix.Transpose(BMatrix);

                // Calculate strain, stress, and internal forces:
                double[] Strain = Elements.Strain(BMatrix, q);
                double[] Stress = Elements.Stress(this.DMatrix, Strain);
                double[] Force  = Elements.InternalForces(BMatrixTranspose, Stress, W, Thickness, detJ);
                TotalInternalForce = VectorAddition(TotalInternalForce, Force);

                // Calculate stiffness:
                double Const = Thickness * W * detJ;
                Matrix Integral = BMatrixTranspose * this.DMatrix;
                Integral *= BMatrix;
                Integral *= Const;
                TotalKMatrix += Integral;
            }
            // Store stiffness and internal forces:
            this.KMatrix = TotalKMatrix;
            this.InternalForce = TotalInternalForce;
        }

        /// <summary>
        /// Calculate Jacobian at (xi,eta) point for 2D Element
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        public static Matrix Jacobian2D(Elements Element, double xi, double eta)
        {
            // Grab nodal locations
            double[] NodalLocations = Element.NodalLocations;

            // Evaluate derivative of shape functions at (xi,eta):
            Matrix dN = ShapeFunctions.EvaluateShapeFunctionDerivatives(Element, xi, eta);

            // Calculate components of Jacobian:
            double[] dX_dXI = dN * NodalLocations;

            // Assemble Jacobian:
            Matrix Jacobian = new(2, 2);
            Jacobian[0, 0] = dX_dXI[0];
            Jacobian[0, 1] = dX_dXI[2];
            Jacobian[1, 0] = dX_dXI[1];
            Jacobian[1, 1] = dX_dXI[3];

            return Jacobian;
        }

        /// <summary>
        /// Calculate B Matrix 2D
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        public static Matrix BMatrix2D(Elements Element, double xi, double eta)
        {
            // Create A matrix (The Stapleton):
            Matrix A = new(3, 4);
            A[0,0] = 1.0; A[1,3] = 1.0; A[2,1] = 1.0; A[2,2] = 1.0;

            // Calculate inverse of jacobian:
            Matrix Jacobian = Elements.Jacobian2D(Element, xi, eta);
            Matrix Jinv = Jacobian.Invert();

            // Assemble inverse jacobian Hat:
            Matrix JInvHat = new(4, 4);
            JInvHat[0,0] = Jinv[0,0]; JInvHat[0,1] = Jinv[0,1];
            JInvHat[1,0] = Jinv[1,0]; JInvHat[1,1] = Jinv[1,1];
            JInvHat[2,2] = Jinv[0,0]; JInvHat[2,3] = Jinv[0,1];
            JInvHat[3,2] = Jinv[1,0]; JInvHat[3,3] = Jinv[1,1];

            // Calculate shape function derivatives
            Matrix dN = ShapeFunctions.EvaluateShapeFunctionDerivatives(Element, xi, eta);

            // Calculate B matrix
            Matrix m = A * JInvHat;
            Matrix BMatrix = m * dN;

            return BMatrix;
        }

        /// <summary>
        /// Calculate strain
        /// </summary>
        /// <param name="B"></B Matrix>
        /// <param name="q"></Nodal Displacements>
        /// <returns></returns>
        public static double[] Strain(Matrix B, double[] q)
        {
            double[] Strain = B * q;
            return Strain;
        }

        /// <summary>
        /// Calculate stress
        /// </summary>
        /// <param name="B"></B Matrix>
        /// <param name="q"></Nodal Displacements>
        /// <returns></returns>
        public static double[] Stress(Matrix D, double[] Strain)
        {
            double[] Stress = D * Strain;
            return Stress;
        }

        /// <summary>
        /// Calculate internal forces
        /// </summary>
        /// <param name="BTranspose"></B matrix transposed>
        /// <param name="Sig"></stress>
        /// <param name="W"></IP weight>
        /// <param name="Thickness"></thickness>
        /// <param name="detJ"></determinant of jacobian>
        /// <returns></returns>
        public static double[] InternalForces(Matrix BTranspose, double[] Sig, double W, double Thickness, double detJ)
        {
            double Const = Thickness * W * detJ;
            double[] InternalForces = BTranspose * Sig;
            InternalForces = InternalForces.Select(r => r * Const).ToArray(); // Multiply internal forces by Const
            return InternalForces;
        }

        /// <summary>
        /// Returns (weights, xi, eta) for each IP
        /// </summary>
        /// <returns></returns>
        public static Matrix GetIntegrationPoints(Elements Element, int NumIPs)
        {

            Matrix W_xi_eta = new(NumIPs, 3); // (weights, xi, eta)
            string type = Element.Type; // element type

            // Square Elements (Quads)
            if (type == "4NQ" || type == "6NQ" || type == "8NQ")
            {
                double W1 = 1.0000000000000000; double W2 = 0.3086419753086425; 
                double W3 = 0.4938271604938276; double W4 = 0.7901234567901236; 
                double A1 = 0.577350269189626;  double A2 = 0.774596669241483;

                switch (NumIPs)
                {
                    case 1:
                        {
                            W_xi_eta[0,0] = 4.0;
                            W_xi_eta[0,1] = 0.0;
                            W_xi_eta[0,2] = 0.0;
                            break;
                        }
                    case 4:
                        {
                            W_xi_eta[0,0] = W1; W_xi_eta[1,0] = W1;  W_xi_eta[2,0] =  W1; W_xi_eta[3,0] = W1;
                            W_xi_eta[0,1] = A1; W_xi_eta[1,1] = A1;  W_xi_eta[2,1] = -A1; W_xi_eta[3,1] = -A1;
                            W_xi_eta[0,2] = A1; W_xi_eta[1,2] = -A1; W_xi_eta[2,2] = A1;  W_xi_eta[3,2] = -A1;
                            break;
                        }
                    case 9:
                        {
                            W_xi_eta[0,0] = W2; W_xi_eta[1,0] = W3;  W_xi_eta[2,0] = W2;  W_xi_eta[3,0] = W3;  W_xi_eta[4,0] = W4;  W_xi_eta[5,0] = W3;  W_xi_eta[6,0] = W2;  W_xi_eta[7,0] = W3;  W_xi_eta[8,0] = W2;
                            W_xi_eta[0,1] = A2; W_xi_eta[1,1] = A2;  W_xi_eta[2,1] = A2;  W_xi_eta[3,1] = 0.0; W_xi_eta[4,1] = 0.0; W_xi_eta[5,1] = 0.0; W_xi_eta[6,1] = -A2; W_xi_eta[7,1] = -A2; W_xi_eta[8,1] = -A2;
                            W_xi_eta[0,2] = A2; W_xi_eta[1,2] = 0.0; W_xi_eta[2,2] = -A2; W_xi_eta[3,2] = A2;  W_xi_eta[4,2] = 0.0; W_xi_eta[5,2] = -A2; W_xi_eta[6,2] = A2;  W_xi_eta[7,2] = 0.0; W_xi_eta[8,2] = -A2;
                            break;
                        }
                    default:
                        {
                            throw new Exception("Incorrect number of IPs for quad element");
                        }
                }
            }
            else if (type == "3NT" || type == "6NT")
            {
                switch (NumIPs)
                {
                    case 1:
                        {
                            W_xi_eta[0,0] = 1.0 / 2.0;
                            W_xi_eta[0,1] = 1.0 / 3.0;
                            W_xi_eta[0,2] = 1.0 / 3.0;
                            break;
                        }
                    case 3:
                        {
                            W_xi_eta[0,0] = 1.0 / 6.0; W_xi_eta[1,0] = 1.0 / 6.0; W_xi_eta[2,0] = 1.0 / 6.0;
                            W_xi_eta[0,1] = 1.0 / 6.0; W_xi_eta[1,1] = 2.0 / 3.0; W_xi_eta[2,1] = 1.0 / 6.0;
                            W_xi_eta[0,2] = 1.0 / 6.0; W_xi_eta[1,2] = 1.0 / 6.0; W_xi_eta[2,2] = 2.0 / 3.0;
                            break;
                        }
                    case 4:
                        {
                            W_xi_eta[0,0] = -9.0 / 32.0; W_xi_eta[1,0] = 25.0 / 96.0; W_xi_eta[2,0] = 25.0 / 96.0; W_xi_eta[3,0] = 25.0 / 96.0;
                            W_xi_eta[0,1] = 1.0 / 3.0;   W_xi_eta[1,1] = 3.0 / 5.0;   W_xi_eta[2,1] = 1.0 / 5.0;   W_xi_eta[3,1] = 1.0 / 5.0;
                            W_xi_eta[0,2] = 1.0 / 3.0;   W_xi_eta[1,2] = 1.0 / 5.0;   W_xi_eta[2,2] = 3.0 / 5.0;   W_xi_eta[3,2] = 1.0 / 5.0;
                            break;
                        }
                    case 7:
                        {
                            W_xi_eta[0,0] = 1.0 / 40.0; W_xi_eta[1,0] = 1.0 / 15.0; W_xi_eta[2,0] = 1.0 / 40.0; W_xi_eta[3,0] = 1.0 / 15.0; W_xi_eta[4,0] = 1.0 / 40.0; W_xi_eta[5,0] = 1.0 / 15.0; W_xi_eta[6,0] = 9.0 / 40.0;
                            W_xi_eta[0,1] = 0.0;        W_xi_eta[1,1] = 1.0 / 2.0;  W_xi_eta[2,1] = 1.0;        W_xi_eta[3,1] = 1.0 / 2.0;  W_xi_eta[4,1] = 0.0;        W_xi_eta[5,1] = 0.0;        W_xi_eta[6,1] = 1.0 / 3.0;
                            W_xi_eta[0,2] = 0.0;        W_xi_eta[1,2] = 0.0;        W_xi_eta[2,2] = 0.0;        W_xi_eta[3,2] = 1.0 / 2.0;  W_xi_eta[4,2] = 1.0;        W_xi_eta[5,2] = 1.0 / 2.0;  W_xi_eta[6,2] = 1.0 / 3.0;
                            break;
                        }
                    default:
                        {
                            throw new Exception("Incorrect number of IPs for triangle element");
                        }
                }
            }
            else if (type == "4NT")
            {
                switch (NumIPs)
                {
                    case 3:
                        {
                            W_xi_eta[0,0] = 5.0 / 36.0;                          W_xi_eta[1,0] = 8.0 / 36.0; W_xi_eta[2,0] = 5.0 / 36.0;
                            W_xi_eta[0,1] = (1.0 / 3.0) - Math.Sqrt(1.0 / 15.0); W_xi_eta[1,1] = 1.0 / 3.0;  W_xi_eta[2,1] = (1.0 / 3.0) + Math.Sqrt(1.0 / 15.0);
                            W_xi_eta[0,2] = 1.0 / 3.0;                           W_xi_eta[1,2] = 1.0 / 3.0;  W_xi_eta[2,2] = 1.0 / 3.0;
                            break;
                        }
                    default:
                        {
                            throw new Exception("4NT element must have 3 IPs");
                        }
                }
            }
            return W_xi_eta;
        }

        /// <summary>
        /// Adds a constant A to a 1D array B
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static double[] VectorAddition(double[] A, double[] B)
        {
            int rowsA = A.Length;
            int rowsB = B.Length;
            if(rowsA != rowsB)
            {
                throw new Exception("Vectors must be same length for addition");
            }

            double[] C = new double[rowsA];
            for (int i = 0; i < rowsA; i++)
            {
                C[i] = A[i] + B[i];
            }
            return C;
        }
    }
}
