using myMath;

namespace FEMAssembly
{
    /// <summary>
    /// Transversly isotropic linear elastic material model
    /// </summary>
    public class TransverslyIsotropicLinearElastic : MaterialModel
    {
        // Constructor
        public TransverslyIsotropicLinearElastic(int NumIPs, int NumStateVars, double E1, double E2, double v23, double G23)
        {
            // State variables
            this.NumStateVars = NumStateVars;
            if (this.NumStateVars > 0)
            {
                StateVars = InitStateVars(this.NumStateVars, NumIPs);
            }
            // Material Properties
            this.E1 = E1;
            this.E2 = E2;
            this.nu23 = v23;
            this.G23 = G23;
        }

        // Methods
        public override void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, int NumIPs, double xi, double eta, double[] Strain, out double[,] DMatrix, out double[] Stress)
        {
            // Calculate DMatrix:
            DMatrix = CalcDMatrixTransverslyIsotropic(E2, nu23, G23);

            // Calculate Stress:
            Stress = Elements.CalcStress(DMatrix, Strain);
        }

        /// <summary>
        /// 2D DMatrix for a transversly isotropic material
        /// </summary>
        public static double[,] CalcDMatrixTransverslyIsotropic(double E2, double nu23, double G23)
        {
            // Compliance matrix
            double[,] SMatrix = new double[3, 3];
            SMatrix[0, 0] = 1.0 / E2;
            SMatrix[0, 1] = -nu23 / E2;
            SMatrix[1, 1] = 1.0 / E2;
            SMatrix[2, 2] = 1.0 / G23;

            // Invert compliance to get stiffness
            double[,] DMatrix = MatrixMath.InvertMatrix(SMatrix);
            return DMatrix;
        }
    }
}
