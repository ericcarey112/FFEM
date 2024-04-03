using System;
using FEMAssembly;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Reflection.Emit;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;

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
        public static double[] SolveDisplacements(Matrix GlobalK, double[] GlobalF)
        {
            Matrix F_Vec = Matrix.VectorToMatrix(GlobalF); // Convert F to Matrix for linear equation solver
            Matrix Q = GlobalK.SolveWith(F_Vec);           // Solve nodal displacements
            double[] GlobalQ = Matrix.MatrixToVector(Q);   // Convert Q from Matrix --> Vector []
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
        public void LoadCheck(Assembly assembly)
        {
            Matrix CurrentLoad = assembly.PeriodicBCs.CurrentMagnitude;
            Matrix TotalLoad = assembly.PeriodicBCs.TotalMagnitude; 
            Matrix LoadDiff = CurrentLoad - TotalLoad;        // Check difference between current load and maximum load
            bool TerminateFlag = Matrix.IsPositive(LoadDiff); // Simulation finishes once CurrentLoad >= MaximumLoad
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
        /// Checks for force equilibrium. If force residual is less than specified tolerance, return converged as true
        /// </summary>
        /// <param name="assembly"></param>
        /// <param name="solver"></param>
        /// <returns></returns>
        public static bool CheckEquilbrium(Assembly assembly, Solver solver)
        {
            // Solve external forces:
            double[] Fext_RVE = Assembly.CalcExternalForces(assembly);

            // Solve internal forces:
            double[] Fint_RVE = new double[assembly.TotalDOF];
            for (int i = 0; i < assembly.NumberOfElements; i++)
            {
                // Add all internal forces from elements:
                Elements element = assembly.ElementList[i];
                double[] Fint_Element = element.InternalForce;
                Fint_RVE = Assembly.LocalToGlobalVector(assembly, element, i, Fint_Element, Fint_RVE);
            }

            // Calculate residual and tolerance:
            solver.Residual = Doubles.SubtractDoubles(Fext_RVE, Fint_RVE);
            double tolerance = Doubles.EuclideanNorm(solver.Residual) / Doubles.EuclideanNorm(Fext_RVE);

            // Check equilibrium:
            bool converged;
            if (tolerance <= solver.ConvergenceTolerance) { converged = true; }
            else { converged = false; }

            return converged;
        }

        /// <summary>
        /// Solves for next displacement guess using newton raphson
        /// </summary>
        public static void ApplyNewtonRaphson(Assembly assembly, double[] Residual)
        {
            // Solve for next displacement guess:
            double[] Q_next = assembly.GlobalK.Invert() * Residual;
            assembly.GlobalQ = Doubles.AddDoubles(assembly.GlobalQ, Q_next);
        }
    }
}
