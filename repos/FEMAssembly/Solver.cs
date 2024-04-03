using myMath;

namespace FEMAssembly
{ 
    /// <summary>
    /// Solver class
    /// </summary>
    public class Solver
    {
        public int SolverType { get; set; }
        public int LoadStepNumber { get; set; }
        public int MaxNRIterations { get; set; }
        public int MaxAttempts { get; set; }
        public int MaxLoadSteps { get; set; }
        public int NRCounter { get; set; }
        public int AttemptCounter { get; set; }
        public bool TerminateFlag { get; set; }
        public double[] Residual {  get; set; }
        public double ConvergenceTolerance { get; set; }

        /// <summary>
        /// Constructor for solver
        /// </summary>
        public Solver()
        {
            this.SolverType = 1;
            this.LoadStepNumber = 0;
            this.MaxNRIterations = 0;
            this.MaxAttempts = 0;
            this.MaxLoadSteps = 0;
            this.NRCounter = 0;
            this.AttemptCounter = 0;
            this.TerminateFlag = false;
            this.Residual = [];
            this.ConvergenceTolerance = 0.0;
        }

        // Methods
        /// <summary>
        /// Solve for nodal dispalcements (Global Q)
        /// </summary>
        /// <param name="GlobalK"></param>
        /// <param name="GlobalF"></param>
        public static double[] SolveDisplacements(double[,] GlobalK, double[] GlobalF)
        {
            double[] GlobalQ = MatrixMath.LinSolve(GlobalK, GlobalF);
            return GlobalQ;
        }

        /// <summary>
        /// Reset attempt counter to 1
        /// </summary>
        public void ResetAttemptCounter()
        {
            this.AttemptCounter = 0;
        }

        /// <summary>
        /// Reset NR counter to 0
        /// </summary>
        public void ResetNRCounter()
        {
            this.NRCounter = 0;
        }

        /// <summary>
        /// Increase attempt counter by 1
        /// </summary>
        public void IncreaseAttemptCounter()
        {
            this.AttemptCounter++;
        }

        /// <summary>
        /// Increase NR Counter
        /// </summary>
        public void IncreaseNRCounter()
        {
            this.NRCounter++;
        }

        /// <summary>
        /// Checks if the max allowable attempts are reached (Number of times load step is cut in half)
        /// </summary>
        public void AttemptCheck()
        {
            if (this.AttemptCounter == this.MaxAttempts)
            {
                this.TerminateFlag = true;
                Console.WriteLine("Maximum number of attempts reached. Simulation ending");
            }
        }

        /// <summary>
        /// Checks current load against max prescribed load
        /// </summary>
        public void LoadCheck(double[,] CurrentLoad, double[,] TotalLoad)
        {
            double[,] LoadDiff = MatrixMath.Subtract(CurrentLoad, TotalLoad); // Check difference between current load and maximum load
            bool TerminateFlag = MatrixMath.IsPositive(LoadDiff); // Simulation finishes once CurrentLoad >= MaximumLoad
            if (TerminateFlag)
            {
                this.TerminateFlag = true;
                Console.WriteLine("Target load reached. Simulation ending");
            }
        }

        /// <summary>
        /// Check if max load steps are reached
        /// </summary>
        public void LoadStepCheck()
        {
            if(this.LoadStepNumber == this.MaxLoadSteps)
            {
                this.TerminateFlag = true;
                Console.WriteLine("Max number of load steps reached. Simulation ending");
            }
        }

        /// <summary>
        /// Solves for next displacement guess using newton raphson (I Want to use Prof Stapleton's NR code, but still don't really understand how the interfaces work... so i'm leaving this as is for now)
        /// </summary>
        public static double[] NewtonRaphson(double[,] GlobalK, double[] OldQ, double[] Residual, string Type)
        {
            double[] NewQ = new double[OldQ.Length];

            if (Type == "Secant")
            {
                // Solve for next displacement guess:
                double[] dQ = MatrixMath.Multiply(MatrixMath.InvertMatrix(GlobalK), Residual);
                NewQ = Doubles.AddDoubles(OldQ, dQ);

            }
            else if (Type == "Tangent")
            {

            }
            else if (Type == "Initial Slope")
            {

            }
            else if (Type == "Numerical")
            {

            }
            else
            {
                throw new Exception("Incorrect Newton-Raphson type in Solver --> NewtonRhapson");
            }

            return NewQ;
        }
    }
}
