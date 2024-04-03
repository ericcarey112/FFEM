using System;
//using FEMAssembly;
using MaterialModels;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEMAssembly
{
    /// <summary>
    /// Structure for fiber element data
    /// </summary>
    public struct FiberData()
    {
        public int ConsitutiveModel { get; set; }
        public double E1 { get; set; }
        public double E2 { get; set; }
        public double v23 { get; set; }
        public double G23 { get; set; }
        public int NumStateVars { get; set; }
    }

    /// <summary>
    /// Structure for matrix element data
    /// </summary>
    public struct MatrixData()
    {
        public int ConsitutiveModel { get; set; }
        public double E2 { get; set; }
        public double v23 { get; set; }
        public double Strength { get; set; }
        public double GIC { get; set; }
        public int NumStateVars { get; set; }
    }

    /// <summary>
    /// Class for reading in input file
    /// </summary>
    public class InputFile
    {
        // Properties
        public string FileName { get; set; }
        public string FileDirectory { get; set; }

        // Constructor
        public InputFile(string filename, string directory)
        {
            this.FileName = filename;
            this.FileDirectory = directory;
        }

        // Methods
        /// <summary>
        /// Generates Input File
        /// </summary>
        public void ReadInputFile(ref Assembly assembly, ref Solver solver, ref ParaviewFiles paraview)
        {
            FiberData fibers = new();
            MatrixData matrix = new();
            PeriodicBCs BCs = new();

            string FilePath = Path.Combine(this.FileDirectory, this.FileName);

            // Make sure file exists
            if (File.Exists(FilePath))
            {
                // Read in file
                using StreamReader reader = new(FilePath);
                string[] lines = File.ReadAllLines(FilePath);

                // Flags
                bool FiberFlag = false;
                bool MatrixFlag = false;
                bool RVEFlag = false;
                bool BCFlag = false;
                bool SolverFlag = false;
                bool ParaviewFlag = false;
                foreach (string line in lines)
                {
                    // Determine current file sections
                    if (string.IsNullOrWhiteSpace(line)) { continue; }
                    if (line == "**FIBERS") { FiberFlag = true; continue; }
                    else if (line == "**MATRIX") { FiberFlag = false; MatrixFlag = true; continue; }
                    else if (line == "**RVE") { MatrixFlag = false; RVEFlag = true; continue; }
                    else if (line == "**BOUNDARY CONDITIONS") { RVEFlag = false; BCFlag = true; continue; }
                    else if (line == "**SOLVER") { BCFlag = false; SolverFlag = true; continue; }
                    else if (line == "**PARAVIEW FILES") { SolverFlag = false; ParaviewFlag = true; continue; }
                    else if (line == "**END") { break; }

                    // Read in data
                    if (FiberFlag) { fibers = ReadFiberData(line, fibers); }
                    if (MatrixFlag) { matrix = ReadMatrixData(line, matrix); }
                    if (RVEFlag) { assembly = ReadRVEData(line, assembly); }
                    if (BCFlag) { BCs = ReadBCData(line, BCs); }
                    if (SolverFlag) { solver = ReadSolverData(line, solver); }
                    if (ParaviewFlag) { paraview = ReadParaview(line, paraview); }
                }

                // Assign material properties and constitutive model to elements:
                for (int i = 0; i < assembly.NumberOfElements; i++)
                {
                    Elements element = assembly.ElementList[i];
                    string type = element.Type;

                    // FIBER ELEMENTS:
                    if (type == "4NT")
                    {
                        // Assign constitutive ModelID to fiber elements:
                        element.ModelID = fibers.ConsitutiveModel;

                        // Isotropic linear elastic constitutive model
                        if (element.ModelID == 1)
                        {
                            element.Model = new IsotropicLinearElastic(element.NumIPs, fibers.E2, fibers.v23);
                        }
                        // Transversly isotropic linear elastic constitutive model
                        else if (element.ModelID == 2)
                        {
                            element.Model = new TransverslyIsotropicLinearElastic(element.NumIPs, fibers.E1, fibers.E2, fibers.v23, fibers.G23);
                        }

                    }

                    // MATRIX ELEMENTS:
                    else if (type == "3NT" || type == "6NT" || type == "6NQ" || type == "8NQ")
                    {
                        // Assign constitutive ModelID to fiber elements:
                        element.ModelID = matrix.ConsitutiveModel;

                        // Isotropic linear elastic constitutive model
                        if (element.ModelID == 1)
                        {
                            element.Model = new IsotropicLinearElastic(element.NumIPs, matrix.E2, matrix.v23);
                        }
                        // Smeared crack constitutive model
                        else if (element.ModelID == 3)
                        {
                            element.Model = new SmearedCrack(element.NumIPs, matrix.E2, matrix.v23, matrix.Strength, matrix.GIC);
                        }
                    }

                    else { throw new Exception("Element type could not be assigned a constitutive model in InputOutput --> ReadInputFile"); }
                }

                // Assign BCs
                assembly.PeriodicBCs = BCs;
                assembly.PeriodicBCs.TotalMagnitude[0, 0] = BCs.E22;
                assembly.PeriodicBCs.TotalMagnitude[0, 1] = BCs.E23;
                assembly.PeriodicBCs.TotalMagnitude[1, 1] = BCs.E33;
            }
        }

        /// <summary>
        /// Read fiber data from input file
        /// </summary>
        private static FiberData ReadFiberData(string line, FiberData fibers)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // Constitutive Model
            if (beforeColon == "ConstitutiveModel")
            {
                if (int.TryParse(afterColon, out int ConstModel))
                {
                    fibers.ConsitutiveModel = ConstModel;
                }
                else
                {
                    throw new Exception("Error in reading fiber constitutive model in InputOutput --> ReadFiberData");
                }
            }

            // E1
            if (beforeColon == "E1")
            {
                if (double.TryParse(afterColon, out double E1))
                {
                    fibers.E1 = E1;
                }
                else
                {
                    throw new Exception("Error in reading fiber E1 data in InputOutput --> ReadFiberData");
                }
            }

            // E2
            if (beforeColon == "E2")
            {
                if (double.TryParse(afterColon, out double E2))
                {
                    fibers.E2 = E2;
                }
                else
                {
                    throw new Exception("Error in reading fiber E2 data in InputOutput --> ReadFiberData");
                }
            }

            // v23
            if (beforeColon == "v23")
            {
                if (double.TryParse(afterColon, out double v23))
                {
                    fibers.v23 = v23;
                }
                else
                {
                    throw new Exception("Error in reading fiber v23 data in InputOutput --> ReadFiberData");
                }
            }

            // G23
            if (beforeColon == "G23")
            {
                if (double.TryParse(afterColon, out double G23))
                {
                    fibers.G23 = G23;
                }
                else
                {
                    throw new Exception("Error in reading fiber G23 data in InputOutput --> ReadFiberData");
                }
            }
            return fibers;
        }

        /// <summary>
        /// Read fiber data from input file
        /// </summary>
        private static MatrixData ReadMatrixData(string line, MatrixData matrix)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // Constitutive Model
            if (beforeColon == "ConstitutiveModel")
            {
                if (int.TryParse(afterColon, out int ConstModel))
                {
                    matrix.ConsitutiveModel = ConstModel;
                }
                else
                {
                    throw new Exception("Error in reading matrix constitutive model in InputOutput --> ReadMatrixData");
                }
            }

            // E2
            if (beforeColon == "E2")
            {
                if (double.TryParse(afterColon, out double E2))
                {
                    matrix.E2 = E2;
                }
                else
                {
                    throw new Exception("Error in reading matrix E2 data in InputOutput --> ReadMatrixData");
                }
            }

            // v23
            if (beforeColon == "v23")
            {
                if (double.TryParse(afterColon, out double v23))
                {
                    matrix.v23 = v23;
                }
                else
                {
                    throw new Exception("Error in reading matrix v23 data in InputOutput --> ReadMatrixData");
                }
            }

            // Strength
            if (beforeColon == "Strength")
            {
                if (double.TryParse(afterColon, out double Strength))
                {
                    matrix.Strength = Strength;
                }
                else
                {
                    throw new Exception("Error in reading matrix strength data in InputOutput --> ReadMatrixData");
                }
            }

            // GIC
            if (beforeColon == "GIC")
            {
                if (double.TryParse(afterColon, out double GIC))
                {
                    matrix.GIC = GIC;
                }
                else
                {
                    throw new Exception("Error in reading matrix GIC data in InputOutput --> ReadMatrixData");
                }
            }
            return matrix;
        }

        /// <summary>
        /// Read in thickness and plane stress or plane strain
        /// </summary>
        private static Assembly ReadRVEData(string line, Assembly assembly)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // Thickness
            if (beforeColon == "Thickness")
            {
                if (double.TryParse(afterColon, out double thickness))
                {
                    assembly.Thickness = thickness;
                }
                else
                {
                    throw new Exception("Error in reading RVE thickness in InputOutput --> ReadRVEData");
                }
            }

            // Plane Stress or Plane Strain
            if (beforeColon == "PlaneStressPlaneStrain")
            {
                if (int.TryParse(afterColon, out int PlaneStressOrStrain))
                {
                    assembly.PlaneStressPlaneStrain = PlaneStressOrStrain;
                }
                else
                {
                    throw new Exception("Error in reading PlaneStressPlaneStrain in InputOutput --> ReadRVEData");
                }
            }
            return assembly;
        }

        /// <summary>
        /// Read in Boundary Conditions
        /// </summary>
        private static PeriodicBCs ReadBCData(string line, PeriodicBCs BCs)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // E22
            if (beforeColon == "E22")
            {
                if (double.TryParse(afterColon, out double E22))
                {
                    BCs.E22 = E22;
                }
                else
                {
                    throw new Exception("Error in reading E22 in InputOutput --> ReadBCData");
                }
            }

            // E33
            if (beforeColon == "E33")
            {
                if (double.TryParse(afterColon, out double E33))
                {
                    BCs.E33 = E33;
                }
                else
                {
                    throw new Exception("Error in reading E33 in InputOutput --> ReadBCData");
                }
            }

            // E23
            if (beforeColon == "E23")
            {
                if (double.TryParse(afterColon, out double E23))
                {
                    BCs.E23 = E23;
                }
                else
                {
                    throw new Exception("Error in reading E23 in InputOutput --> ReadBCData");
                }
            }

            // Load Step
            if (beforeColon == "LoadStep")
            {
                if (double.TryParse(afterColon, out double LoadStep))
                {
                    Matrix InitLoad = new(2, 2);
                    if (BCs.E22 != 0.0) { InitLoad[0, 0] = LoadStep; }
                    if (BCs.E23 != 0.0) { InitLoad[0, 1] = LoadStep; }
                    if (BCs.E33 != 0.0) { InitLoad[1, 1] = LoadStep; }
                    BCs.InitialLoadStep = InitLoad;
                }
                else
                {
                    throw new Exception("Error in reading LoadStep in InputOutput --> ReadBCData");
                }
            }
            return BCs;
        }

        /// <summary>
        /// Read in solver settings
        /// </summary>
        private static Solver ReadSolverData(string line, Solver solver)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // Solver Type
            if (beforeColon == "Type")
            {
                if (int.TryParse(afterColon, out int type))
                {
                    solver.SolverType = type;
                }
                else
                {
                    throw new Exception("Error in reading solver type in InputOutput --> ReadSolverData");
                }
            }

            // Max NR Iterations
            if (beforeColon == "NRMax")
            {
                if (int.TryParse(afterColon, out int NRMax))
                {
                    solver.MaxNRIterations = NRMax;
                }
                else
                {
                    throw new Exception("Error in reading NRMax in InputOutput --> ReadSolverData");
                }
            }

            // Max Attempts
            if (beforeColon == "MaxAttempts")
            {
                if (int.TryParse(afterColon, out int MaxAttempts))
                {
                    solver.MaxAttempts = MaxAttempts;
                }
                else
                {
                    throw new Exception("Error in reading MaxAttempts in InputOutput --> ReadSolverData");
                }
            }

            // Max number of load steps
            if (beforeColon == "MaxLoadSteps")
            {
                if (int.TryParse(afterColon, out int MaxLoadSteps))
                {
                    solver.MaxLoadSteps = MaxLoadSteps;
                }
                else
                {
                    throw new Exception("Error in reading MaxLoadSteps in InputOutput --> ReadSolverData");
                }
            }

            // Convergence criteria
            if (beforeColon == "ConvDelta")
            {
                if (double.TryParse(afterColon, out double ConvDelta))
                {
                    solver.ConvergenceTolerance = ConvDelta;
                }
                else
                {
                    throw new Exception("Error in reading ConvDelta in InputOutput --> ReadSolverData");
                }
            }
            return solver;
        }

        /// <summary>
        /// Read in whether or not to generate paraview vtks
        /// </summary>
        private static ParaviewFiles ReadParaview(string line, ParaviewFiles paraview)
        {
            string[] parts = line.Split(':');
            string beforeColon = parts[0].Trim();
            string afterColon = parts[1].Trim();

            // Solver Type
            if (beforeColon == "ParaviewVTKs")
            {
                if (int.TryParse(afterColon, out int GenerateVTKs))
                {
                    if (GenerateVTKs == 1) { paraview.GenerateFiles = true; }
                    else { paraview.GenerateFiles = false; }
                }
                else
                {
                    throw new Exception("Error in reading ParaviewVTKs in InputOutput --> ReadParaview");
                }
            }
            return paraview;
        }
    }
}
