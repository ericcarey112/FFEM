using myMath;

namespace FEMAssembly
{
    /// <summary>
    /// Abstract 2D Element Class
    /// </summary>
    public abstract class Elements
    {
        // Properties
        public int ModelID { get; set; }
        public int ElementNumber { get; set; }
        public string Type { get; set; }
        public int NumDim { get; set; }
        public int NumNodes { get; set; }
        public int NDOFPNode { get; set; }
        public int TotalDOF {  get; set; }
        public int NumIPs {  get; set; }
        public MaterialModel Model { get; set; }
        public double Thickness { get; set; }
        public int PlaneStressPlaneStrain { get; set; }
        public double[] NodalLocations { get; set; }
        public double[] NodalDisplacements { get; set; }
        public double[] InternalForce { get; set; }
        public double[] ForceVector { get; set; }
        public double[,] KMatrix { get; set; }

        // Constructor (How do i format a constructor for an abstract class? Do i make a method to initialize these for each inherited element class?)
        public Elements()
        {
            this.ModelID = 1;
            this.ElementNumber = 0;
            this.Type = " ";
            this.NumDim = 0;
            this.NumNodes = 0;
            this.NDOFPNode = 0;
            this.TotalDOF = 0;
            this.NumIPs = 1;
            this.Model = new IsotropicLinearElastic(this.NumIPs, 0, 0.0, 0.0);
            this.Thickness = 1.0;
            this.PlaneStressPlaneStrain = 2;
            this.NodalLocations = new double[1];
            this.NodalDisplacements = new double[1];
            this.InternalForce = new double[1];
            this.ForceVector = new double[1];
            this.KMatrix = new double[1, 1];
        }

        /// <summary>
        /// Calculate Jacobian at (xi,eta) point for 2D Element
        /// </summary>
        /// <param name="Element"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        public static double[,] CalcJacobian2D(string type, double[] NodalLocations, double xi, double eta)
        {
            // Evaluate derivative of shape functions at (xi,eta):
            double[,] dN = ShapeFunctions.EvaluateShapeFunctionDerivatives(type, xi, eta);

            // Calculate components of Jacobian:
            double[] dX_dXI = MatrixMath.Multiply(dN, NodalLocations);

            // Assemble Jacobian:
            double[,] Jacobian = new double[2, 2];
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
        public static double[,] BMatrix2D(string type, double[] NodalLocations, double xi, double eta)
        {
            // Create A matrix (The Stapleton):
            double[,] A = new double[3, 4];
            A[0,0] = 1.0; A[1,3] = 1.0; A[2,1] = 1.0; A[2,2] = 1.0;

            // Calculate inverse of jacobian:
            double[,] Jacobian = CalcJacobian2D(type, NodalLocations, xi, eta);
            double[,] Jinv = MatrixMath.InvertMatrix(Jacobian);

            // Assemble inverse jacobian Hat:
            double[,] JInvHat = new double[4, 4];
            JInvHat[0,0] = Jinv[0,0]; JInvHat[0,1] = Jinv[0,1];
            JInvHat[1,0] = Jinv[1,0]; JInvHat[1,1] = Jinv[1,1];
            JInvHat[2,2] = Jinv[0,0]; JInvHat[2,3] = Jinv[0,1];
            JInvHat[3,2] = Jinv[1,0]; JInvHat[3,3] = Jinv[1,1];

            // Calculate shape function derivatives
            double[,] dN = ShapeFunctions.EvaluateShapeFunctionDerivatives(type, xi, eta);

            // Calculate B matrix
            double[,] m = MatrixMath.Multiply(A, JInvHat);
            double[,] BMatrix = MatrixMath.Multiply(m, dN);

            return BMatrix;
        }

        /// <summary>
        /// Calculate strain
        /// </summary>
        /// <param name="B"></B Matrix>
        /// <param name="q"></Nodal Displacements>
        /// <returns></returns>
        public static double[] CalcStrain(double[,] B, double[] q)
        {
            double[] Strain = MatrixMath.Multiply(B, q);
            return Strain;
        }

        /// <summary>
        /// Calculate stress
        /// </summary>
        /// <param name="D"></D Matrix>
        /// <param name="Strain"></Strain>
        /// <returns></returns>
        public static double[] CalcStress(double[,] DMatrix, double[] Strain)
        {
            double[] Stress = MatrixMath.Multiply(DMatrix, Strain);
            return Stress;
        }

        /// <summary>
        /// Calculate principal stress or strain
        /// </summary>
        /// <param name=""></param>
        /// <param name="theata"></param>
        /// <returns></returns>
        public static double[] RotateStressOrStrain(double[] A, double theta)
        {
            double[] RotatedStressOrStrain = new double[A.Length];

            // Components of vector
            double Axx = A[0];
            double Ayy = A[1];
            double Axy = A[2];

            // Assembly A into a matrix:
            double[,] A_Tens = new double[3, 3];
            A_Tens[0, 0] = Axx; A_Tens[0, 1] = Axy;
            A_Tens[1, 0] = Axy; A_Tens[1, 1] = Ayy;
            A_Tens[2, 2] = 1.0;

            // Assemble rotation matrix and its transpose:
            double[,] R = new double[3, 3];
            R[0, 0] =  Math.Cos(theta); R[0, 1] = Math.Sin(theta);
            R[1, 0] = -Math.Sin(theta); R[1, 1] = Math.Cos(theta);
            R[2, 2] = 1.0;
            
            double[,] RTranspose = MatrixMath.Transpose(R);

            // Calculate rotated components:
            double[,] M1 = MatrixMath.Multiply(A_Tens, R);
            double[,] PrincTens = MatrixMath.Multiply(RTranspose, M1);
            RotatedStressOrStrain[0] = PrincTens[0, 0];
            RotatedStressOrStrain[1] = PrincTens[1, 1];
            RotatedStressOrStrain[2] = 0.0;

            return RotatedStressOrStrain;
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
        public static double[] CalcInternalForces(double[,] BTranspose, double[] Stress, double W, double Thickness, double detJ)
        {
            double Const = Thickness * W * detJ;
            double[] InternalForces = MatrixMath.Multiply(BTranspose, Stress);
            InternalForces = InternalForces.Select(r => r * Const).ToArray(); // Multiply internal forces by Const
            return InternalForces;
        }

        /// <summary>
        /// Function to call for calculating characteristic length around an IP
        /// </summary>
        public static double CalcCharLength(string type, double[] NodalLocations, int NumIPs, double theta, double xi, double eta)
        {
            double CharLength;
            // Triangular elements
            if (type == "3NT" || type == "4NT" || type == "6NT")
            {
                CharLength = CalcCharLength_Tri(type, NodalLocations, theta, xi, eta);
            }

            // Quadrilateral elements
            else if (type == "4NQ" || type == "6NQ" || type == "8NQ")
            {
                CharLength = CalcCharLength_Quad(type, NodalLocations, NumIPs, theta, xi, eta);
            }

            else
            {
                throw new Exception("Incorrect Element Type for Calculating Charactersitic Length Elements --> CalcCharLength");
            }

            return CharLength;
        }

        /// <summary>
        /// Calculates the characteristic length of a 2D quad element
        /// </summary>
        private static double CalcCharLength_Quad(string type, double[] NodalLocations, int NumIPs, double theta, double xi, double eta)
        {
            // Map theta from global to isoparametric space:
            theta = MapThetaFromGlobalToIso(type, NodalLocations, theta, xi, eta);

            // Create grid lines around integration points:
            double[] Lines_xi;
            double[] Lines_eta;
            switch (NumIPs)
            {
                case 1:
                    {
                        Lines_xi  = [-1.0, 1.0];
                        Lines_eta = [-1.0, 1.0];
                        break;
                    }
                case 4:
                    {
                        Lines_xi  = [-1.0, 0.0, 1.0];
                        Lines_eta = [-1.0, 0.0, 1.0];
                        break;
                    }
                case 9:
                    {
                        Lines_xi  = [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0];
                        Lines_eta = [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0];
                        break;
                    }
                default:
                    {
                        throw new Exception("Incorrect number of IPs Elements -> CalcCharLength_Quad");
                    }
            }

            // Calculate distances from IP to each line:
            int count = 0;
            double[] dists = new double[Lines_xi.Length + Lines_eta.Length];
            double[] xipts = new double[dists.Length];
            double[] etapts = new double[dists.Length];
            for (int i = 0; i < Lines_xi.Length; i++)
            {
                // Find (xi,eta) along grid line:
                xipts[count] = Lines_xi[i]; 
                etapts[count] = -1.0 / Math.Tan(theta) * (xipts[count] - xi) + eta;

                // Calculate distance between (xi,eta) of IP and (xi eta) along grid line:
                double[] p1 = [xi, eta];
                double[] p2 = [xipts[count], etapts[count]];
                dists[count] = CalcDistBetweenPoints(p1, p2);
                count++;
            }

            for (int i = 0; i < Lines_eta.Length; i++)
            {
                // Find (xi,eta) along grid line:
                etapts[count] = Lines_eta[i];
                xipts[count] = -(etapts[count] - eta) * Math.Tan(theta) + xi;

                // Calculate distance between (xi,eta) of IP and (xi,eta) along grid line:
                double[] p1 = [xi, eta];
                double[] p2 = [xipts[count], etapts[count]];
                dists[count] = CalcDistBetweenPoints(p1, p2);
                count++;
            }

            // Find minumum two distances and corresponding (xi,eta) points:
            int minIndex1 = 0;
            int minIndex2 = 0;
            for (int i = 0; i < dists.Length; i++)
            {
                if (dists[i] < dists[minIndex1]) {  minIndex1 = i; }
            }
            for (int i = 0; i < dists.Length; i++)
            {
                if (i != minIndex1 && (dists[i] < dists[minIndex2]) || minIndex1 == minIndex2) { minIndex2 = i; }
            }
            List<double[]> Points_xieta = [[xipts[minIndex1], etapts[minIndex1]], [xipts[minIndex2], etapts[minIndex2]]];

            // Map points from (xi,eta) to (x,y)
            double[,] XY = new double[2, 2];
            for (int i = 0; i < 2; i++)
            {
                double[] N = ShapeFunctions.EvaluateShapeFunctions(type, Points_xieta[i][0], Points_xieta[i][1]);
                double[,] NMat = ShapeFunctions.AssembleShapeFunctions(type, N);
                double[] xy = MatrixMath.Multiply(NMat, NodalLocations);
                XY[i, 0] = xy[0];
                XY[i, 1] = xy[1];
            }

            // Calculate distance between points (This is the characteristic length)
            double[] xy1 = [XY[0, 0], XY[0, 1]];
            double[] xy2 = [XY[1, 0], XY[1, 1]];
            double CharLength = CalcDistBetweenPoints(xy1, xy2);
            return CharLength;
        }

        /// <summary>
        /// Calculates the characteristic length of a 2D triangular element
        /// </summary>
        private static double CalcCharLength_Tri(string type, double[] NodalLocations, double theta, double xi, double eta)
        {
            // Map theta from global to isoparametric space:
            theta = MapThetaFromGlobalToIso(type, NodalLocations, theta, xi, eta);

            // Create lines around integration point:
            // Hypotenuse
            double xi_hyp = (xi / Math.Tan(theta) + eta - 1.0) * Math.Pow((1.0 / Math.Tan(theta) - 1.0), -1);
            double eta_hyp = -xi_hyp + 1.0;
            // Height
            double xi_height = 0.0;
            double eta_height = -1.0 / Math.Tan(theta) * (0.0 - xi) + eta;
            // Base
            double xi_base = -(0.0 - eta) * Math.Tan(theta) + xi;
            double eta_base = 0.0;
            // Assemble points:
            List<double[]> IntPoints = [[xi_hyp, eta_hyp], [xi_height, eta_height], [xi_base, eta_base]];

            // Calculate distance between IP and lines around IP:
            double[] p1 = [xi, eta];
            double[] dists = new double[3];
            double[] p2;
            for (int i = 0; i < 3; i++)
            {
                p2 = [IntPoints[i][0], IntPoints[i][1]];
                dists[i] = CalcDistBetweenPoints(p1, p2);
            }

            // Find minumum two distances and corresponding (xi,eta) points:
            int minIndex1 = 0;
            int minIndex2 = 0;
            for (int i = 0; i < dists.Length; i++)
            {
                if (dists[i] < dists[minIndex1]) { minIndex1 = i; }
            }
            for (int i = 0; i < dists.Length; i++)
            {
                if (i != minIndex1 && (dists[i] < dists[minIndex2]) || minIndex1 == minIndex2) { minIndex2 = i; }
            }
            List<double[]> Points_xieta =
            [
                [IntPoints[minIndex1][0], IntPoints[minIndex1][1]],
                [IntPoints[minIndex2][0], IntPoints[minIndex2][1]],
            ];

            // Map points from (xi,eta) to (x,y)
            double[,] XY = new double[2, 2];
            for (int i = 0; i < 2; i++)
            {
                double[] N = ShapeFunctions.EvaluateShapeFunctions(type, Points_xieta[i][0], Points_xieta[i][1]);
                double[,] NMat = ShapeFunctions.AssembleShapeFunctions(type, N);
                double[] xy = MatrixMath.Multiply(NMat, NodalLocations);
                XY[i, 0] = xy[0];
                XY[i, 1] = xy[1];
            }

            // Calculate distance between points (This is the characteristic length)
            double[] xy1 = [XY[0, 0], XY[0, 1]];
            double[] xy2 = [XY[1, 0], XY[1, 1]];
            double CharLength = CalcDistBetweenPoints(xy1, xy2);
            return CharLength;
        }

        /// <summary>
        /// Maps angle theta from global to isoparametric space
        /// </summary>
        private static double MapThetaFromGlobalToIso(string type, double[] NodalLocations, double theta, double xi, double eta)
        {
            // Components of vector w.r.t theta
            double[] xy = [1.0, Math.Tan(theta)];

            // Calculate Jacobian at (xi,eta):
            double[,] Jac = CalcJacobian2D(type, NodalLocations, xi, eta);

            // Calculate theta in isoparametric:
            double[] xieta = MatrixMath.Multiply(Jac, xy);
            xi = xieta[0];
            eta = xieta[1];
            theta = Math.Atan(eta / xi);

            return theta;
        }

        private static double CalcDistBetweenPoints(double[] p1, double[] p2)
        {
            double dx = p2[0] - p1[0];
            double dy = p2[1] - p1[1];
            double distance = Math.Sqrt(dx * dx + dy * dy);
            return distance;
        }

        /// <summary>
        /// Returns (weights, xi, eta) for each IP
        /// </summary>
        /// <returns></returns>
        public static double[,] GetIntegrationPoints(Elements Element)
        {
            int NumIPs = Element.NumIPs; // Number of IPs
            string type = Element.Type;  // Element type

            // Square Elements (Quads)
            double[,] W_xi_eta = new double[NumIPs, 3]; // (weights, xi, eta)
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
    }
}
