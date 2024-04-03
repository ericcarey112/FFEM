namespace FEMAssembly
{
    /// <summary>
    /// Abstract class for creating a material model
    /// </summary>
    public abstract class MaterialModel
    {
        // Properties
        public int NumStateVars;
        public double[,] StateVars = new double[0, 0];
        public double E1;
        public double E2;
        public double nu12;
        public double nu23;
        public double G23;

        /// <summary>
        /// Solves for stiffness and internal forces based on the constitutive model
        /// </summary>
        public abstract void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, int NumIPs, double xi, double eta, double[] Strain, out double[,] DMatrix, out double[] Stress);

        /// <summary>
        /// Initialize state variables
        /// </summary>
        /// <param name="NumStateVars"></param>
        /// <param name="NumIPs"></param>
        public static double[,] InitStateVars(int NumStateVars, int NumIPs)
        {
            double[,] StateVars = new double[NumStateVars, NumIPs];
            for (int i = 0; i < NumStateVars; i++)
            {
                for (int j = 0; j < NumIPs; j++)
                {
                    StateVars[i, j] = 0.0;
                }
            }
            return StateVars;
        }
    }
}
