using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MaterialModels
{
    /// <summary>
    /// Smeared crack material model
    /// </summary>
    public class SmearedCrack : MaterialModel
    {
        // Properties
        public double Strength { get; set; }
        public double GIC { get; set; }

        // Constructor
        public SmearedCrack(int NumIPs, double YoungMod, double PoissonRatio, double Strength, double GIC)
        {
            // State variables
            this.NumStateVars = 3;
            if (this.NumStateVars > 0)
            {
                InitStateVars(this.NumStateVars, NumIPs);
            }
            // Material properties
            this.E2 = YoungMod;
            this.nu23 = PoissonRatio;
            this.Strength = Strength;
            this.GIC = GIC;
        }

        // Methods
        public override void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, double xi, double eta, double[] Strain, out Matrix DMatrix, out double[] Stress)
        {
            // Set max allowable damage:
            double DMax = 0.9999;

            // Get state variables:
            double NormToCrackAngle = StateVars[IPNum][0]; // Angle normal to crack
            double CurrentDamage = StateVars[IPNum][1];    // Damage
            double InitStrain = StateVars[IPNum][2];       // Initiation Strain

            // Calculate current strength:
            double CurrentStrength = (1.0 - CurrentDamage) * Strength;

            // Calculate DMatrix:
            DMatrix = CalcDMatrixIsotropicDamage(E2, nu23, PlaneStressPlaneStrain, CurrentDamage);

            // Calculate stress:
            Stress = CalcStress(DMatrix, Strain);

            // Calculate normal to crack angle:
            if (CurrentDamage == 0.0)
            {
                StateVars[IPNum][0] = CalcNormToCrackAngle(Stress);
            }

            // Assemble transformation matrix and reuter matrix:
            Matrix T = AssembleTransformationMatrix(NormToCrackAngle);
            Matrix R = AssembleReuterMatrix();
            Matrix Tinv = T.Invert();
            Matrix Rinv = R.Invert();

            // Rotate stresses, strains, and stiffness to the principal frame:
            double[] RotatedStrains = Elements.RotateStressOrStrain(Strain, NormToCrackAngle);
            double[] RotatedStresses = Elements.RotateStressOrStrain(Stress, NormToCrackAngle);
            RotateDMatrixGlobalToPrincipal(ref DMatrix, T, Tinv, R, Rinv);

            // Determine max norm to crack stress and corresponding strain:
            double NormToCrackStress = 0.0;
            double NormToCrackStrain = 0.0;
            int direction = 0;
            for (int i = 0; i < Strain.Length; i++)
            {
                if (RotatedStresses[i] > NormToCrackStress)
                {
                    NormToCrackStress = RotatedStresses[i];
                    NormToCrackStrain = RotatedStrains[i];
                    direction = i;
                }
            }

            // Check direction of max stress normal to crack surface:
            if (CurrentDamage == 0)
            {
                if (direction == 2) { NormToCrackAngle += Math.PI / 2.0; }
            }

            // Calculate characteristic length:
            //double CharLength = Elements.CalcCharLength(element, NormToCrackAngle, xi, eta);
            double CharLength = 1.0;

            // Calculate failure strain:
            double FailStrain = 2.0 * GIC / (Strength * CharLength);

            // Check for damage initiation:
            if (NormToCrackStress >= CurrentStrength && CurrentDamage == 0.0)
            {
                // Calculate initiation strain:
                InitStrain = NormToCrackStrain * Strength / NormToCrackStress;
                StateVars[IPNum][0] = NormToCrackAngle;
                StateVars[IPNum][2] = InitStrain;
            }

            // Check for damage progression:
            if (NormToCrackStress >= CurrentStrength && CurrentDamage < DMax)
            {
                // Calculate damage:
                double UpdatedDamage = 1.0 - (((1.0 - FailStrain / NormToCrackStrain) * InitStrain) / (InitStrain - FailStrain));

                // Make sure UpdatedDamage < DMax
                if (UpdatedDamage >= DMax) { UpdatedDamage = DMax; }

                // Damage can't be negative
                if (UpdatedDamage < 0.0) { UpdatedDamage = 0.0; }

                // Store damage:
                StateVars[IPNum][1] = UpdatedDamage;

                // Recalculate DMatrix with new damage and rotate back to global frame:
                DMatrix = CalcDMatrixIsotropicDamage(E2, nu23, element.PlaneStressPlaneStrain, CurrentDamage);
                RotateDMatrixPrincipalToGlobal(ref DMatrix, T, Tinv, R, Rinv);

                // Recalculate stress:
                Stress = Elements.CalcStress(DMatrix, Strain);
            }
        }

        /// <summary>
        /// Rotate the DMatrix from global to principal frame
        /// </summary>
        /// <returns></returns>
        private static void RotateDMatrixGlobalToPrincipal(ref Matrix DMatrix, Matrix T, Matrix Tinv, Matrix R, Matrix Rinv)
        {
            Matrix M1 = T * DMatrix;
            Matrix M2 = M1 * R;
            Matrix M3 = M2 * Tinv;
            DMatrix = M3 * Rinv;
        }

        private static void RotateDMatrixPrincipalToGlobal(ref Matrix DMatrix, Matrix T, Matrix Tinv, Matrix R, Matrix Rinv)
        {
            Matrix M1 = Tinv * DMatrix;
            Matrix M2 = M1 * R;
            Matrix M3 = M2 * T;
            DMatrix = M3 * Rinv;
        }

        /// <summary>
        /// Calculate angle normal to crack
        /// </summary>
        /// <param name="model"></param>
        private static double CalcNormToCrackAngle(double[] Stress)
        {
            // Components of stress
            double Sxx = Stress[0];
            double Syy = Stress[1];
            double Sxy = Stress[2];

            double theta;
            // Pure shear case (45 degrees exactly):
            if (Math.Abs(Sxx - Syy) < 0.0000001)
            {
                theta = Math.PI / 4.0;
            }
            // All other cases:
            else
            {
                theta = Math.Atan(2.0 * Sxy / (Sxx - Syy)) / 2.0;
            }

            return theta;
        }

        /// <summary>
        /// Assembles reuter matrix for converting between voigt and tensorial shear strains
        /// </summary>
        /// <returns></returns>
        private static Matrix AssembleReuterMatrix()
        {
            Matrix R = new(3, 3);
            R[0, 0] = 1.0; R[0, 1] = 0.0; R[0, 2] = 0.0;
            R[1, 0] = 0.0; R[1, 1] = 1.0; R[1, 2] = 0.0;
            R[2, 0] = 0.0; R[2, 1] = 0.0; R[2, 2] = 2.0;
            return R;
        }

        /// <summary>
        /// Transformation matrix (Two rotation matrices multiplied together)
        /// </summary>
        /// <param name="theta"></param>
        /// <returns></returns>
        private static Matrix AssembleTransformationMatrix(double theta)
        {
            Matrix Trans = new(3, 3);
            Trans[0, 0] = Math.Pow(Math.Cos(theta), 2); Trans[0, 1] = Math.Pow(Math.Sin(theta), 2); Trans[0, 2] = 2.0 * Math.Sin(theta) * Math.Cos(theta);
            Trans[1, 0] = Math.Pow(Math.Sin(theta), 2); Trans[1, 1] = Math.Pow(Math.Cos(theta), 2); Trans[1, 2] = -2.0 * Math.Sin(theta) * Math.Cos(theta);
            Trans[2, 0] = Math.Sin(theta) * Math.Cos(theta); Trans[2, 1] = -Math.Sin(theta) * Math.Cos(theta); Trans[2, 2] = Math.Pow(Math.Cos(theta), 2) - Math.Pow(Math.Sin(theta), 2);
            return Trans;
        }

        /// <summary>
        /// 2D DMatrix for isotropic plane stress or strain element with isotropic damage
        /// </summary>
        public static Matrix CalcDMatrixIsotropicDamage(double E, double nu, int PlaneStressPlaneStrain, double damage)
        {
            Matrix DMatrix = new(3, 3);

            // Calculate D (Plane Stress):
            if (PlaneStressPlaneStrain == 1)
            {
                double constant = E / (1.0 - nu * nu);
                DMatrix[0, 0] = (1.0 - damage) * constant;
                DMatrix[0, 1] = (1.0 - damage) * nu * constant;
                DMatrix[0, 2] = 0.0;
                DMatrix[1, 0] = DMatrix[0, 1];
                DMatrix[1, 1] = DMatrix[0, 0];
                DMatrix[1, 2] = 0.0;
                DMatrix[2, 0] = DMatrix[0, 2];
                DMatrix[2, 1] = DMatrix[1, 2];
                DMatrix[2, 2] = (1.0 - damage) * constant * (1 - nu) / (2.0);
            }
            // Calculate D (Plane Strain):
            else if (PlaneStressPlaneStrain == 2)
            {
                double constant = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                DMatrix[0, 0] = (1.0 - damage) * constant * (1.0 - nu);
                DMatrix[0, 1] = (1.0 - damage) * constant * nu;
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

        /// <summary>
        /// 2D DMatrix for isotropic plane stress or strain element with isotropic damage
        /// </summary>
        public static Matrix CalcDMatrixAnisotropicDamage(double E, double nu, int PlaneStressPlaneStrain, double damage)
        {

            Matrix DMatrix = new(3, 3);

            // Calculate D (Plane Stress):
            if (PlaneStressPlaneStrain == 1)
            {
                double constant = E / (1.0 - nu * nu);
                DMatrix[0, 0] = (1.0 - damage) * constant;
                DMatrix[0, 1] = nu * constant;
                DMatrix[0, 2] = 0.0;
                DMatrix[1, 0] = DMatrix[0, 1];
                DMatrix[1, 1] = constant;
                DMatrix[1, 2] = 0.0;
                DMatrix[2, 0] = DMatrix[0, 2];
                DMatrix[2, 1] = DMatrix[1, 2];
                DMatrix[2, 2] = (1.0 - damage) * constant * (1 - nu) / (2.0);
            }
            // Calculate D (Plane Strain):
            else if (PlaneStressPlaneStrain == 2)
            {
                double constant = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                DMatrix[0, 0] = (1.0 - damage) * constant * (1.0 - nu);
                DMatrix[0, 1] = constant * nu;
                DMatrix[0, 2] = 0.0;
                DMatrix[1, 0] = DMatrix[0, 1];
                DMatrix[1, 1] = constant * (1.0 - nu);
                DMatrix[1, 2] = 0.0;
                DMatrix[2, 0] = DMatrix[0, 2];
                DMatrix[2, 1] = DMatrix[1, 2];
                DMatrix[2, 2] = constant * (1 - 2.0 * nu) / (2.0);
            }
            return DMatrix;
        }
    }
}
