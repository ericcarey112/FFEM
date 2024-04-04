using myMath;

namespace FEMAssembly
{ 
    /// <summary>
    /// Dirichlet (Displacements applied directly to nodes)
    /// </summary>
    public struct DirichletBCs()
    {
        public int NumberOfDimensions = 2;
        public int Direction = 0;
        public double InitialLoadStep = 0.0;
        public double CurrentLoadStep = 0.0;
        public double CurrentMagnitude = 0.0;
        public double TotalMagnitude = 0.0;
    }

    /// <summary>
    /// Periodic (Periodic displacments applied from localizing far field strains)
    /// </summary>
    public struct PeriodicBCs()
    {
        public int NumberOfDimensions = 2;
        public double E22 = 0.0;
        public double E23 = 0.0;
        public double E33 = 0.0;
        public double[,] InitialLoadStep = new double[2, 2];
        public double[,] CurrentLoadStep = new double[2, 2];
        public double[,] CurrentMagnitude = new double[2, 2];
        public double[,] TotalMagnitude = new double[2, 2];
    }

    public class Assembly // Should this be an asbtract class with different analysis types? Like implicit/explicit inherited classes for example...) Should ask Stapleton about this
    {
        // Properties
        public int NumberOfElements { get; set; }
        public int NumberOfNodes { get; set; }
        public int NDOFPNode { get; set; }
        public int TotalDOF { get; set; }
        public List<int[]> NodePairArray { get; set; }
        public List<double[]> NodalLocations { get; set; }
        public List<int[]> Connectivity {  get; set; }
        public double[,] GlobalK { get; set; }
        public double[] GlobalF { get; set; }
        public double[] GlobalQ { get; set; }
        public double[,] GlobalM { get; set; }
        public double[,] K_PreBC { get; set; }
        public double[] F_PreBC { get; set; }
        public double[,] K_PreBC_Init { get; set; }
        public double[] F_PreBC_Init { get; set; }
        public double PenaltyStiffness { get; set; }
        public DirichletBCs DirichletBCs { get; set; }
        public PeriodicBCs PeriodicBCs { get; set; }
        public List<Elements> ElementList { get; set; }
        public List<int> RightEdgeNodes { get; set; }
        public List<int> TopEdgeNodes { get; set; }
        public int PinnedNode { get; set; }
        public int PlaneStressPlaneStrain {  get; set; }
        public double Thickness { get; set; }
        public double[] RVELength { get; set; }
        public double[] HomogenizedStress { get; set; }
        public double[] HomogenizedStrain { get; set; }
        public double[,] HomogenizedStiffness { get; set; }
        public double[] ExternalForces { get; set; }
        public double[] InternalForces {  get; set; }
        private double[] Residual { get; set; }

        // Constructor
        public Assembly() // Most of these properties (NumberOfElements, NumberofNodes, etc get assigned during the mesh reading. Kind of awkward to have a bunch of "... = 0" for a bunch of properties. Best way to set up constructors?) Should ask Stapleton
        {
            this.NumberOfElements = 0;
            this.NumberOfNodes = 0;
            this.NDOFPNode = 2;
            this.TotalDOF = 0;
            this.NodePairArray = [];
            this.NodalLocations = [];
            this.Connectivity = new List<int[]>();
            this.GlobalK = new double[1, 1];
            this.GlobalF = new double[1];
            this.GlobalQ = new double[1];
            this.GlobalM = new double[1, 1];
            this.K_PreBC = new double[1, 1];
            this.F_PreBC = new double[1];
            this.K_PreBC_Init = new double[1, 1];
            this.F_PreBC_Init = new double [1];
            this.PenaltyStiffness = 0.0;
            this.DirichletBCs = new DirichletBCs();
            this.PeriodicBCs  = new PeriodicBCs();
            this.ElementList = [];
            this.RightEdgeNodes = [];
            this.TopEdgeNodes = [];
            this.PinnedNode = 0;
            this.PlaneStressPlaneStrain = 2;
            this.Thickness = 0.0;
            this.RVELength = new double[2];
            this.HomogenizedStress = new double [3];
            this.HomogenizedStrain = new double [3];
            this.HomogenizedStiffness = new double[3, 3];
            this.ExternalForces = [];
            this.InternalForces = [];
            this.Residual = [];
        }

        // Methods
        /// <summary>
        /// Submit model to solver
        /// </summary>
        public void AnalyzeRVE(Solver solver, OutputFile output, ParaviewFiles paraview)
        {
            // Loop over each load step:
            for (int i = 0; i < solver.MaxLoadSteps; i++)
            {
                // Check if model should terminate:
                if (solver.TerminateFlag) { break; }

                // Increase load step:
                IncreaseLoadStep();
                solver.LoadStepNumber++;

                // Solve constitutive laws and look for force equilibrium:
                bool converged = false;
                while (converged == false)
                {
                    // Loop over each element and apply constitutive model:
                    for (int j = 0; j < NumberOfElements; j++)
                    {
                        // Preallocate stiffness and internal forces for each element:
                        Elements element = ElementList[j];
                        element.KMatrix = new double[element.TotalDOF, element.TotalDOF];
                        element.InternalForce = new double[element.TotalDOF];

                        // Assign global nodal displacements to each element:
                        element.NodalDisplacements = GlobalToLocalVector(element, Connectivity, j, GlobalQ);

                        // Get integration points for each element:
                        double[,] W_xi_eta = Elements.GetIntegrationPoints(element);

                        // Calculate stiffness and internal force at each IP:
                        for (int k = 0; k < element.NumIPs; k++)
                        {
                            // Get Weights and IPs:
                            double W   = W_xi_eta[k, 0];
                            double xi  = W_xi_eta[k, 1];
                            double eta = W_xi_eta[k, 2];

                            // Calculate Jacobian and determinant:
                            double[,] Jacobian = Elements.CalcJacobian2D(element.Type, element.NodalLocations, xi, eta);
                            double detJ = MatrixMath.Determinant(Jacobian);

                            // Calculate B Matrix and its transpose:
                            double[,] BMatrix = Elements.BMatrix2D(element.Type, element.NodalLocations, xi, eta);
                            double[,] BMatrixTranspose = MatrixMath.Transpose(BMatrix);

                            // Calculate strain:
                            double[] Strain = Elements.CalcStrain(BMatrix, element.NodalDisplacements);

                            // Calculate DMatrix and stress with material model:
                            element.Model.SolveDMatrixAndStress(element.Type, element.NodalLocations, element.PlaneStressPlaneStrain, k, element.NumIPs, xi, eta, Strain, out double[,] DMatrix, out double[] Stress);

                            // Calculate stiffness at each IP:
                            double Const = element.Thickness * W * detJ;
                            double[,] Integral = MatrixMath.Multiply(BMatrixTranspose, DMatrix);
                            Integral = MatrixMath.Multiply(Integral, BMatrix);
                            Integral = MatrixMath.ScalarMultiply(Const, Integral);
                            element.KMatrix = MatrixMath.Add(element.KMatrix, Integral);

                            // Calculate internal force at each IP:
                            double[] Force = Elements.CalcInternalForces(BMatrixTranspose, Stress, W, element.Thickness, detJ);
                            element.InternalForce = Doubles.AddDoubles(element.InternalForce, Force);
                        }
                    }

                    // Assemble global K and F:
                    AssembleGlobalK();
                    AssembleGlobalF();

                    // Store initial K and F at beginning of load step and apply BCs:
                    if (solver.NRCounter == 0 && solver.AttemptCounter == 0)
                    {
                        K_PreBC_Init = K_PreBC;
                        F_PreBC_Init = F_PreBC;
                        ApplyPeriodicBCs();
                    }

                    // Check convergence:
                    if (solver.SolverType == 1)      // Static linear solver
                    {
                        converged = true;
                    }
                    else if (solver.SolverType == 2) // Static non-linear solver
                    {
                        // Calculate Residual (Fext - Fint)
                        CalcExternalForcesRVE();
                        AssembleInternalForcesRVE();
                        Residual = Doubles.SubtractDoubles(ExternalForces, InternalForces);

                        // Check Equilibrium:
                        converged = CheckEquilibrium(ExternalForces, Residual, solver.ConvergenceTolerance);
                        
                    }
                    else { throw new Exception("Invalid solver type. Must be 1 (Static linear) or 2 (Static non-linear)"); }

                    // If converged, post process and move to next load step:
                    if (converged)
                    {
                        // Assemble force sum matrix on first load step:
                        if (solver.LoadStepNumber == 1) { AssembleGlobalM(); }

                        // Calculate homogenized stress, strain, and stiffness:
                        CalcHomogenizedStressStrain();
                        CalcHomogenizedStiffness();

                        // Generate output files:
                        output.UpdateOutputFile(this, solver); // Homogenzied data
                        if (paraview.GenerateFiles) { paraview.GenerateVTKFiles(this, solver.LoadStepNumber); } // Paraview files

                        // Display converged to command line:
                        Console.WriteLine("Load Step " + solver.LoadStepNumber + " Converged. Attempts: " + solver.AttemptCounter + " NR Iterations: " + solver.NRCounter + " Max Load: " + MatrixMath.GetMax(PeriodicBCs.CurrentMagnitude));

                        // Solver checks
                        solver.LoadCheck(PeriodicBCs.CurrentMagnitude, PeriodicBCs.TotalMagnitude); // Specified load reached? (Total Magnitude)
                        solver.LoadStepCheck(); // Max number of load steps reached? 

                        // Reset NR and Attempt Counter:
                        solver.ResetNRCounter();
                        solver.ResetAttemptCounter();
                        continue;
                    }

                    // If not converged, apply Newton-Raphson:
                    else
                    {
                        // If NR Counter = MaxNRIterations, cut load step in half:
                        if (solver.NRCounter == solver.MaxNRIterations)
                        {
                            // Check if max attempts are reached:
                            solver.IncreaseAttemptCounter();
                            solver.AttemptCheck();
                            if (solver.TerminateFlag) { return; }

                            // Reset NR Counter:
                            solver.ResetNRCounter();

                            // Reset K and F to beginning of load step:
                            K_PreBC = K_PreBC_Init;
                            F_PreBC = F_PreBC_Init;

                            // Reduce load step
                            ReduceLoadStep();
                        }
                        // Use NR to guess next displacements:
                        GlobalQ = Solver.NewtonRaphson(GlobalK, GlobalQ, Residual, "Secant");
                        solver.IncreaseNRCounter();
                    }
                }
            }
        }

        /// <summary>
        /// Adds all the element objects into a list (ONLY FOR SIMPLE 8NQ UNIT TESTS --- NOT FOR ACTUAL CODE)
        /// </summary>
        /// <param name="NumberOfElements"></param>
        /// <param name="NodalLocations"></param>
        public void AddElementsToList(int NumberOfElements, List<double[]> NodalLocations)
        {
            for (int i = 0; i < NumberOfElements; i++)
            {
                this.ElementList.Add(new Element_8NQ());
                this.ElementList[i].ElementNumber = i + 1;
                this.ElementList[i].NodalLocations = NodalLocations[i];
                this.ElementList[i].PlaneStressPlaneStrain = this.PlaneStressPlaneStrain;
            }
            
        }
        
        /// <summary>
        /// Assemble Global K
        /// </summary>
        private void AssembleGlobalK()
        {
            double[,] GlobalK = new double[TotalDOF, TotalDOF];

            for (int i = 1; i < NumberOfElements+1; i++)
            {
                Elements element = ElementList[i - 1]; // Current element
                for(int j = 1; j < element.NumNodes+1; j++)
                {
                    int CurrentNode = (int)this.Connectivity[i-1][j-1]; // Current node number
                    for(int k = 1; k < element.NumNodes+1; k++)
                    {
                        for (int l = 1; l < element.NumDim+1; l++)
                        {
                            for (int m = 1; m < element.NumDim+1; m++)
                            {
                                // Global K
                                GlobalK[(element.NumDim * CurrentNode - l), (element.NumDim * (int)this.Connectivity[i - 1][k - 1] - m)] = GlobalK[(element.NumDim * CurrentNode - l), (element.NumDim * (int)this.Connectivity[i - 1][k - 1] - m)] + element.KMatrix[(element.NumDim * j - l), (element.NumDim * k - m)];
                            }
                        }
                    }
                }
            }
            K_PreBC = GlobalK;
        }

        /// <summary>
        /// Assemble Global F
        /// </summary>
        private void AssembleGlobalF()
        {
            double[] GlobalF = new double[TotalDOF];

            for (int i = 1; i < NumberOfElements + 1; i++)
            {
                Elements element = ElementList[i - 1]; // Current element
                for (int j = 1; j < element.NumNodes + 1; j++)
                {
                    int CurrentNode = (int)this.Connectivity[i - 1][j - 1]; // Current node number
                    for (int k = 1; k < element.NumDim + 1; k++)
                    {
                        // Global F
                        GlobalF[element.NumDim * CurrentNode - k] = GlobalF[element.NumDim * CurrentNode - k] + element.ForceVector[element.NumDim * j - k];
                    }
                }
            }
            F_PreBC = GlobalF;
        }

        /// <summary>
        /// Assembles force sum matrix (M) for homogenization
        /// </summary>
        private void AssembleGlobalM()
        {
            double[,] GlobalM = new double[NDOFPNode * NDOFPNode, TotalDOF];

            // Tractions on right face (Txx and Txy):
            for (int i = 0; i < RightEdgeNodes.Count; i++)
            {
                for(int j = 0; j < NDOFPNode; j++)
                {
                    int DOF = NDOFPNode * RightEdgeNodes[i] + (j - NDOFPNode);
                    GlobalM[j, DOF] = 1;
                } 
            }

            // Tractions on top face (Tyx and Tyy):
            for (int i = 0; i < TopEdgeNodes.Count; i++)
            {
                for (int j = 0; j < NDOFPNode; j++)
                {
                    int DOF = NDOFPNode * TopEdgeNodes[i] + (j - NDOFPNode);
                    GlobalM[NDOFPNode * NDOFPNode + (j - NDOFPNode) , DOF] = 1;
                }
            }
            this.GlobalM = GlobalM;
        }

        /// <summary>
        /// Apply periodic BCs
        /// </summary>
        private void ApplyPeriodicBCs()
        {
            // Inputs
            double[,] FarFieldStrain = PeriodicBCs.CurrentMagnitude;
            double[,] GlobalK = K_PreBC;
            double[] GlobalF = F_PreBC;

            // Apply BCs using penalty method:
            double PenaltyStiffness = MatrixMath.GetMax(GlobalK) * Math.Pow(10, 6);
            for (int i = 0; i < NodePairArray.Count; i++)
            {
                // Get nodes:
                int n1 = NodePairArray[i][0];
                int n2 = NodePairArray[i][1];

                // Get displacements from far field strain:
                double n1x = NodalLocations[n1-1][0];
                double n1y = NodalLocations[n1-1][1];
                double n2x = NodalLocations[n2-1][0];
                double n2y = NodalLocations[n2-1][1];
                double[] ub1 = GetBoundaryDispalcements(n1x, n1y, FarFieldStrain);
                double[] ub2 = GetBoundaryDispalcements(n2x, n2y, FarFieldStrain);

                for(int j = 0; j < NDOFPNode; j++)
                {
                    // Get global index of the dof and nodes:
                    int I1 = NDOFPNode * n1 + (j - NDOFPNode);
                    int I2 = NDOFPNode * n2 + (j - NDOFPNode);

                    // Apply to appropriate entries of F:
                    GlobalF[I1] = GlobalF[I1] + PenaltyStiffness * (ub1[j] - ub2[j]);
                    GlobalF[I2] = GlobalF[I2] + PenaltyStiffness * (ub2[j] - ub1[j]);

                    // Apply to appropriate entries of K:
                    GlobalK[I1, I1] = GlobalK[I1, I1] + PenaltyStiffness;
                    GlobalK[I2, I2] = GlobalK[I2, I2] + PenaltyStiffness;
                    GlobalK[I1, I2] = GlobalK[I1, I2] - PenaltyStiffness;
                    GlobalK[I2, I1] = GlobalK[I2, I1] - PenaltyStiffness;
                }
            }

            // Apply penalty method to pinned node:
            for(int i = 0; i < NDOFPNode; i++)
            {
                int I1 = NDOFPNode * PinnedNode + (i - NDOFPNode);
                GlobalK[I1, I1] = GlobalK[I1, I1] + PenaltyStiffness;
            }
            this.GlobalF = GlobalF;
            this.GlobalK = GlobalK;
            this.PenaltyStiffness = PenaltyStiffness;
        }

        /// <summary>
        /// Calcuate boundary displacements from far field strains
        /// </summary>
        private static double[] GetBoundaryDispalcements(double x, double y, double[,] FarFieldStrain)
        {
            // Inputs
            double[,] H = new double[2, 2];
            double[] xy = [x, y];

            H[0, 0] = FarFieldStrain[0, 0];
            H[0, 1] = 2.0 * FarFieldStrain[0, 1];
            H[1, 0] = FarFieldStrain[1, 0];
            H[1, 1] = FarFieldStrain[1, 1];

            double[] ub = MatrixMath.Multiply(H, xy);

            return ub;
        }

        /// <summary>
        /// Calculate external forces of RVE
        /// </summary>
        private void CalcExternalForcesRVE()
        {
            double[] Fext = new double[TotalDOF];

            // Inputs
            double[,] FarFieldStrain = PeriodicBCs.CurrentMagnitude;
  
            for (int i = 0; i < NodePairArray.Count; i++)
            {
                // Get the nodes:
                int n1 = NodePairArray[i][0];
                int n2 = NodePairArray[i][1];

                // Get displacements from far field strain:
                double n1x = NodalLocations[n1 - 1][0];
                double n1y = NodalLocations[n1 - 1][1];
                double n2x = NodalLocations[n2 - 1][0];
                double n2y = NodalLocations[n2 - 1][1];
                double[] ub1 = GetBoundaryDispalcements(n1x, n1y, FarFieldStrain);
                double[] ub2 = GetBoundaryDispalcements(n2x, n2y, FarFieldStrain);

                // Determine if shear loading (Multiply displacements by 2)
                if (FarFieldStrain[0, 1] > 0.0)
                {
                    ub1 = Doubles.MultiplyVectorByConstant(ub1, 2.0);
                    ub2 = Doubles.MultiplyVectorByConstant(ub2, 2.0);
                }

                for (int j = 0; j < NDOFPNode; j++)
                {
                    // DOFs
                    int DOF_N1 = n1 * NDOFPNode + (j - NDOFPNode);
                    int DOF_N2 = n2 * NDOFPNode + (j - NDOFPNode);

                    // Calculate external forces
                    Fext[DOF_N1] += -PenaltyStiffness * (GlobalQ[DOF_N1] - GlobalQ[DOF_N2]) + PenaltyStiffness * (ub1[j] - ub2[j]);
                    Fext[DOF_N2] += -PenaltyStiffness * (GlobalQ[DOF_N2] - GlobalQ[DOF_N1]) + PenaltyStiffness * (ub2[j] - ub1[j]);
                }
            }
            this.ExternalForces = Fext;
        }

        /// <summary>
        /// Assemble the internal forces of RVE
        /// </summary>
        private void AssembleInternalForcesRVE()
        {
            double[] Fint_RVE = new double[TotalDOF];

            for (int j = 0; j < NumberOfElements; j++)
            {
                // Add all internal forces from elements:
                Elements element = ElementList[j];
                Fint_RVE = LocalToGlobalVector(element, Connectivity, j, element.InternalForce, Fint_RVE);
            }
            this.InternalForces = Fint_RVE;
        }

        /// <summary>
        /// Check if force residual is less than specified tolerance
        /// </summary>
        /// <returns></returns>
        private static bool CheckEquilibrium(double[] force, double[] Residual, double minTol)
        {
            // Calculate tolerance:
            //double tol = Doubles.EuclideanNorm(Residual) / Doubles.EuclideanNorm(force);
            double tol = Math.Sqrt(Math.Abs(Doubles.MultiplyDoubles(Residual, Residual)));
 
            // Check calculated tolerance (tol) vs. minimum specified tolerance (minTol)
            if (tol <= minTol) { return true; }
            else { return false; }
        }

        /// <summary>
        /// Calculate homogenized stress and strain of an RVE
        /// </summary>
        public void CalcHomogenizedStressStrain()
        {
            // Inputs
            double[,] K_PreBC = this.K_PreBC;
            double[,] GlobalM = this.GlobalM;
            double[] GlobalQ = this.GlobalQ;
            double thickness = this.Thickness;
            double Ly = this.RVELength[0];
            double Lz = this.RVELength[1];

            // Assemble inverse of areas matrix:
            double[,] Ainv = CalcInvAreas(thickness, Ly, Lz);

            // Calculate homogenized stress:
            double[] A1 = MatrixMath.Multiply(K_PreBC, GlobalQ);
            double[] A2 = MatrixMath.Multiply(GlobalM, A1);
            this.HomogenizedStress = MatrixMath.Multiply(Ainv, A2);

            // Calculate homogenized strain:
            this.HomogenizedStrain[0] = this.PeriodicBCs.CurrentMagnitude[0, 0];
            this.HomogenizedStrain[1] = this.PeriodicBCs.CurrentMagnitude[1, 1];
            this.HomogenizedStrain[2] = this.PeriodicBCs.CurrentMagnitude[0, 1];
        }

        /// <summary>
        /// Calcualates all homogenized stiffness components (See slides 6_6_2023 for full derivation)
        /// </summary>
        public void CalcHomogenizedStiffness()
        {
            double[,] F_E = new double[TotalDOF, 3];

            // Calculate derivatives w.r.t. strain
            for (int i = 0;  i < NodePairArray.Count; i++)
            {
                // Get the nodes:
                int n1 = NodePairArray[i][0];
                int n2 = NodePairArray[i][1];

                // Get the nodal locations:
                double n1x = NodalLocations[n1 - 1][0];
                double n1y = NodalLocations[n1 - 1][1];
                double n2x = NodalLocations[n2 - 1][0];
                double n2y = NodalLocations[n2 - 1][1];

                for(int j = 0; j < NDOFPNode; j++)
                {
                    // DOFs
                    int DOF_N1 = NDOFPNode * n1 + (j - NDOFPNode);
                    int DOF_N2 = NDOFPNode * n2 + (j - NDOFPNode);

                    // Derivative of force w.r.t. strain:
                    int mod = j % 2;
                    if (mod == 0) // index j is even (Recall indexing starts at 0. This is why corresponding matlab is if mod != 0)
                    {
                        F_E[DOF_N1, 0] += PenaltyStiffness * (n1x - n2x);
                        F_E[DOF_N1, 2] += PenaltyStiffness * (n1y - n2y);
                        F_E[DOF_N2, 0] += PenaltyStiffness * (n2x - n1x);
                        F_E[DOF_N2, 2] += PenaltyStiffness * (n2y - n1y);
                    }
                    else // index j is odd
                    {
                        F_E[DOF_N1, 1] += PenaltyStiffness * (n1y - n2y);
                        F_E[DOF_N2, 1] += PenaltyStiffness * (n2y - n1y);
                    }
                }
            }

            // Calculate homogenized stiffness:
            double[,] Ainv = CalcInvAreas(Thickness, RVELength[0], RVELength[1]);
            double[,] K_Inv = MatrixMath.InvertMatrix(GlobalK);
            double[,] A1 = MatrixMath.Multiply(K_Inv, F_E);
            double[,] A2 = MatrixMath.Multiply(K_PreBC, A1);
            double[,] A3 = MatrixMath.Multiply(GlobalM, A2);
            this.HomogenizedStiffness = MatrixMath.Multiply(Ainv, A3);
        }

        /// <summary>
        /// Inverse areas matrix used for stress and stiffness homogenization
        /// </summary>
        /// <param name="Thickness"></param>
        /// <param name="Ly"></param>
        /// <param name="Lz"></param>
        /// <returns></returns>
        private static double[,] CalcInvAreas(double Thickness, double Ly, double Lz)
        {
            double[,] InvAreas = new double[3, 4];
            double Axy = Thickness * Ly;
            double Axz = Thickness * Lz;

            InvAreas[0, 0] = 1.0 / Axy;
            InvAreas[1, 3] = 1.0 / Axz;
            InvAreas[2, 1] = 1.0 / (2.0 * Axy);
            InvAreas[2, 2] = 1.0 / (2.0 * Axz);

            return InvAreas;
        }

        /// <summary>
        /// Increases the load step by the prescribed initial load step or a reduced load step that is cut in half
        /// </summary>
        /// <param name="solver"></param>
        public void IncreaseLoadStep()
        {
            PeriodicBCs ModifiedBCs = this.PeriodicBCs;

            ModifiedBCs.CurrentLoadStep = ModifiedBCs.InitialLoadStep;
            ModifiedBCs.CurrentMagnitude = MatrixMath.Add(ModifiedBCs.CurrentMagnitude, ModifiedBCs.CurrentLoadStep);
            this.PeriodicBCs = ModifiedBCs;
        }

        /// <summary>
        /// Reduces the load step by half if NRCounter = MaxNRIterations
        /// </summary>
        public void ReduceLoadStep()
        {
            PeriodicBCs ModifiedBCs = this.PeriodicBCs;

            // Reset load to beginning of load step and increase by half:
            ModifiedBCs.CurrentMagnitude = MatrixMath.Subtract(ModifiedBCs.CurrentMagnitude, ModifiedBCs.CurrentLoadStep);
            ModifiedBCs.CurrentLoadStep = MatrixMath.ScalarMultiply(0.5, ModifiedBCs.CurrentLoadStep);
            ModifiedBCs.CurrentMagnitude = MatrixMath.Add(ModifiedBCs.CurrentMagnitude, ModifiedBCs.CurrentLoadStep);
            this.PeriodicBCs = ModifiedBCs;
        }

        /// <summary>
        /// Distributes a global vector (Assembly level) to a local vector (Element level) 
        /// </summary>
        /// <returns></returns>
        public static double[] GlobalToLocalVector(Elements element, List<int[]> ConnectMat, int index, double[] GlobalVec)
        {
            double[] LocalVec = new double[element.TotalDOF];

            // Assign local vector based on connectivity matrix:
            for (int i = 0; i < element.NumNodes; i++)
            {
                int Node = ConnectMat[index][i];
                for (int j = 0; j < element.NDOFPNode; j++)
                {
                    int DOF_Global = element.NumDim * Node - (element.NumDim - j);
                    int DOF_Local = element.NumDim * i + j;                
                    LocalVec[DOF_Local] = GlobalVec[DOF_Global];
                }

            }
            return LocalVec;
        }

        /// <summary>
        /// Adds a local vector (Element level) to a global vector (Assembly level)
        /// </summary>
        /// <returns></returns>
        public static double[] LocalToGlobalVector(Elements element, List<int[]> ConnectMat, int index, double[] LocalVec, double[] GlobalVec)
        {
            // Add to global vector based on connectivity matrix:
            for (int i = 0; i < element.NumNodes; i++)
            {
                int Node = ConnectMat[index][i];
                for(int j = 0; j < element.NDOFPNode; j++)
                {
                    int DOF_Global = element.NumDim * Node - (element.NumDim - j);
                    int DOF_Local = element.NumDim * i + j;
                    GlobalVec[DOF_Global] = GlobalVec[DOF_Global] + LocalVec[DOF_Local];
                }
            }
            return GlobalVec;
        }
    }
}
