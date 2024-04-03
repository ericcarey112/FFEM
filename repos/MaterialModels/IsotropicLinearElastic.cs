using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MaterialModels
{
    /// <summary>
    /// Isotropic linear elastic material model
    /// </summary>
    public class IsotropicLinearElastic : MaterialModel
    {
        // Constructor
        public IsotropicLinearElastic(int NumIPs, double YoungMod, double PoissonRatio)
        {
            // State variables
            this.NumStateVars = 0;
            if (this.NumStateVars > 0)
            {
                InitStateVars(this.NumStateVars, NumIPs);
            }
            // Material properties
            this.E2 = YoungMod;
            this.nu23 = PoissonRatio;
        }

        // Methods
        public override void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, double xi, double eta, double[] Strain, out Matrix DMatrix, out double[] Stress)
        {
            // Calculate DMatrix
            DMatrix = CalcDMatrixIsotropic(E2, nu23, PlaneStressPlaneStrain);

            // Calculate stress:
            Stress = CalcStress(DMatrix, Strain);
        }

        // Methods
        /// <summary>
        /// 2D DMatrix for isotropic plane stress or strain material
        /// </summary>
        public static Matrix CalcDMatrixIsotropic(double E, double nu, int PlaneStressPlaneStrain)
        {
            Matrix DMatrix = new(3, 3);

            // Calculate D (Plane Stress):
            if (PlaneStressPlaneStrain == 1)
            {
                double constant = E / (1.0 - nu * nu);
                DMatrix[0, 0] = constant;
                DMatrix[0, 1] = nu * constant;
                DMatrix[0, 2] = 0.0;
                DMatrix[1, 0] = DMatrix[0, 1];
                DMatrix[1, 1] = DMatrix[0, 0];
                DMatrix[1, 2] = 0.0;
                DMatrix[2, 0] = DMatrix[0, 2];
                DMatrix[2, 1] = DMatrix[1, 2];
                DMatrix[2, 2] = constant * (1 - nu) / (2.0);
            }
            // Calculate D (Plane Strain):
            else if (PlaneStressPlaneStrain == 2)
            {
                double constant = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                DMatrix[0, 0] = constant * (1.0 - nu);
                DMatrix[0, 1] = constant * nu;
                DMatrix[0, 2] = 0.0;
                DMatrix[1, 0] = DMatrix[0, 1];
                DMatrix[1, 1] = DMatrix[0, 0];
                DMatrix[1, 2] = 0.0;
                DMatrix[2, 0] = DMatrix[0, 2];
                DMatrix[2, 1] = DMatrix[1, 2];
                DMatrix[2, 2] = constant * (1 - 2.0 * nu) / (2.0);
            }
            return DMatrix;
        }
    }
}
