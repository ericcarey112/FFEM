using System.ComponentModel;

namespace FEMAssembly
{
    /// <summary>
    /// Abstract class for creating a material model
    /// </summary>
    public abstract class MaterialModel
    {
        // Properties
        public int NumStateVars;
        public List<double[]> StateVars  = [];
        public double E1;
        public double E2;
        public double nu12;
        public double nu23;
        public double G23;

        /// <summary>
        /// Solves for stiffness and internal forces based on the constitutive model
        /// </summary>
        public abstract void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, double xi, double eta, double[] Strain, out Matrix DMatrix, out double[] Stress);

        /// <summary>
        /// Initialize state variables
        /// </summary>
        /// <param name="NumStateVars"></param>
        /// <param name="NumIPs"></param>
        public void InitStateVars(int NumStateVars, int NumIPs)
        {
            double[] zero = new double[NumStateVars];
            for (int i = 0; i < NumIPs; i++)
            {
                StateVars.Add(zero);
            }
        }
        /// <summary>
        /// Calculate stress for any material model
        /// </summary>
        /// <param name="DMatrix"></param>
        /// <param name="Strain"></param>
        /// <returns></returns>
        public static double[] CalcStress(Matrix DMatrix, double[] Strain)
        {
            double[] Stress;
            return Stress = DMatrix * Strain;
        }
    }
}
