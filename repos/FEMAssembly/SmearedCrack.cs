using myMath;

namespace FEMAssembly
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
        public SmearedCrack(int NumIPs, int NumStateVars, double YoungMod, double PoissonRatio, double Strength, double GIC)
        {
            // State variables
            this.NumStateVars = NumStateVars;
            if (this.NumStateVars > 0)
            {
                StateVars = InitStateVars(this.NumStateVars, NumIPs);
            }
            // Material properties
            this.E2 = YoungMod;
            this.nu23 = PoissonRatio;
            this.Strength = Strength;
            this.GIC = GIC;
        }

        // Methods
        public override void SolveDMatrixAndStress(string type, double[] NodalLocations, int PlaneStressPlaneStrain, int IPNum, int NumIPs, double xi, double eta, double[] Strain, out double[,] DMatrix, out double[] Stress)
        {
            // Set max allowable damage:
            double DMax = 0.9999;

            // Get state variables:
            double NormToCrackAngle = StateVars[0, IPNum]; // Angle normal to crack
            double CurrentDamage = StateVars[1, IPNum];    // Damage
            double InitStrain = StateVars[2, IPNum];       // Initiation Strain

            // Calculate current strength:
            double CurrentStrength = (1.0 - CurrentDamage) * Strength;

            // Calculate DMatrix:
            DMatrix = CalcDMatrixIsotropicDamage(E2, nu23, PlaneStressPlaneStrain, CurrentDamage);

            // Calculate stress:
            Stress = Elements.CalcStress(DMatrix, Strain);

            // Calculate normal to crack angle:
            if (CurrentDamage == 0.0)
            {
                NormToCrackAngle = CalcNormToCrackAngle(Stress);
                StateVars[0, IPNum] = NormToCrackAngle;
            }

            // Assemble transformation matrix and reuter matrix:
            double[,] T = AssembleTransformationMatrix(NormToCrackAngle);
            double[,] R = AssembleReuterMatrix();
            double[,] Tinv = MatrixMath.InvertMatrix(T);
            double[,] Rinv = MatrixMath.InvertMatrix(R);

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
            double CharLength = Elements.CalcCharLength(type, NodalLocations, NumIPs, NormToCrackAngle, xi, eta);

            // Calculate failure strain:
            double FailStrain = 2.0 * GIC / (Strength * CharLength);

            // Check for damage initiation:
            if (NormToCrackStress >= CurrentStrength && CurrentDamage == 0.0)
            {
                // Calculate initiation strain:
                InitStrain = NormToCrackStrain * Strength / NormToCrackStress;
                StateVars[0, IPNum] = NormToCrackAngle;
                StateVars[2, IPNum] = InitStrain;
                
            }

            // Check for damage progression:
            if (NormToCrackStress >= CurrentStrength && CurrentDamage < DMax)
            {
                // Calculate damage:
                double UpdatedDamage = 1.0 - (((1.0 - FailStrain / NormToCrackStrain) * InitStrain) / (InitStrain - FailStrain));

                // Damage can't decrease:
                if (UpdatedDamage < CurrentDamage) { UpdatedDamage = CurrentDamage; }

                // Make sure UpdatedDamage < DMax
                if (UpdatedDamage >= DMax) { UpdatedDamage = DMax; }

                // Damage can't be negative
                if (UpdatedDamage < 0.0) { UpdatedDamage = 0.0; }

                // Store damage:
                StateVars[1, IPNum] = UpdatedDamage;

                // Recalculate DMatrix with new damage and rotate back to global frame:
                DMatrix = CalcDMatrixIsotropicDamage(E2, nu23, PlaneStressPlaneStrain, UpdatedDamage);
                RotateDMatrixPrincipalToGlobal(ref DMatrix, T, Tinv, R, Rinv);

                // Recalculate stress:
                Stress = Elements.CalcStress(DMatrix, Strain);
            }
        }

        /// <summary>
        /// Rotate the DMatrix from global to principal frame
        /// </summary>
        /// <returns></returns>
        private static void RotateDMatrixGlobalToPrincipal(ref double[,] DMatrix, double[,] T, double[,] Tinv, double[,] R, double[,] Rinv)
        {
            double[,] M1 = MatrixMath.Multiply(T, DMatrix);
            double[,] M2 = MatrixMath.Multiply(M1, R);
            double[,] M3 = MatrixMath.Multiply(M2, Tinv);
            DMatrix = MatrixMath.Multiply(M3, Rinv);
        }

        private static void RotateDMatrixPrincipalToGlobal(ref double[,] DMatrix, double[,] T, double[,] Tinv, double[,] R, double[,] Rinv)
        {
            double[,] M1 = MatrixMath.Multiply(Tinv, DMatrix);
            double[,] M2 = MatrixMath.Multiply(M1, R);
            double[,] M3 = MatrixMath.Multiply(M2, T);
            DMatrix = MatrixMath.Multiply(M3, Rinv);
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
        private static double[,] AssembleReuterMatrix()
        {
            double[,] R = new double[3, 3];
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
        private static double[,] AssembleTransformationMatrix(double theta)
        {
            double[,] Trans = new double[3, 3];
            Trans[0, 0] = Math.Pow(Math.Cos(theta), 2); Trans[0, 1] = Math.Pow(Math.Sin(theta), 2); Trans[0, 2] = 2.0 * Math.Sin(theta) * Math.Cos(theta);
            Trans[1, 0] = Math.Pow(Math.Sin(theta), 2); Trans[1, 1] = Math.Pow(Math.Cos(theta), 2); Trans[1, 2] = -2.0 * Math.Sin(theta) * Math.Cos(theta);
            Trans[2, 0] = Math.Sin(theta) * Math.Cos(theta); Trans[2, 1] = -Math.Sin(theta) * Math.Cos(theta); Trans[2, 2] = Math.Pow(Math.Cos(theta), 2) - Math.Pow(Math.Sin(theta), 2);
            return Trans;
        }

        /// <summary>
        /// 2D DMatrix for isotropic plane stress or strain element with isotropic damage
        /// </summary>
        public static double[,] CalcDMatrixIsotropicDamage(double E, double nu, int PlaneStressPlaneStrain, double damage)
        {
            double[,] DMatrix = new double[3, 3];

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
        public static double[,] CalcDMatrixAnisotropicDamage(double E, double nu, int PlaneStressPlaneStrain, double damage)
        {

            double[,] DMatrix = new double[3, 3];

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
