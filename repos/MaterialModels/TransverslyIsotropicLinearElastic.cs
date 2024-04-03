using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MaterialModels
{
    /// <summary>
    /// Transversly isotropic linear elastic material model
    /// </summary>
    public class TransverslyIsotropicLinearElastic : MaterialModel
    {
        // Constructor
        public TransverslyIsotropicLinearElastic(int NumIPs, double E1, double E2, double v23, double G23)
        {
            // State variables
            this.NumStateVars = 0;
            if (this.NumStateVars > 0)
            {
                InitStateVars(this.NumStateVars, NumIPs);
            }
            // Material Properties
            this.E1 = E1;
            this.E2 = E2;
            this.nu23 = v23;
            this.G23 = G23;
        }

        // Methods
        public override void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, double xi, double eta, double[] Strain, out Matrix DMatrix, out double[] Stress)
        {
            // Calculate DMatrix:
            DMatrix = CalcDMatrixTransverslyIsotropic(E2, nu23, G23);

            // Calculate Stress:
            Stress = CalcStress(DMatrix, Strain);
        }

        /// <summary>
        /// 2D DMatrix for a transversly isotropic material
        /// </summary>
        public static Matrix CalcDMatrixTransverslyIsotropic(double E2, double nu23, double G23)
        {
            // Compliance matrix
            Matrix SMatrix = new(3, 3);
            SMatrix[0, 0] = 1.0 / E2;
            SMatrix[0, 1] = -nu23 / E2;
            SMatrix[1, 1] = 1.0 / E2;
            SMatrix[2, 2] = 1.0 / G23;

            // Invert compliance to get stiffness
            Matrix DMatrix = SMatrix.Invert();
            return DMatrix;
        }
    }
}
