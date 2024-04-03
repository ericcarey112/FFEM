namespace FEMAssembly
{
    class Program
    {
        /// <summary>
        /// Main program
        /// </summary>
        /// <param name="args"></param>
        public static void Main()
        {
            ///////////////// INPUTS /////////////////
            // Create instance of each class:
            MeshFile mesh = new("Hex_20Pack.txt", "C:\\Users\\Eric_Carey.UMLADCO\\OneDrive - UMass Lowell\\FFEM_CSharp\\MeshFiles\\RVEs");           // (Filename.txt, Directory)
            InputFile inputfile = new("ExampleInputFile.txt", "C:\\Users\\Eric_Carey.UMLADCO\\OneDrive - UMass Lowell\\FFEM_CSharp\\InputFiles");    // (Filename.txt, Directory)
            OutputFile outputfile = new("Hex_20Pack.txt", "C:\\Users\\Eric_Carey.UMLADCO\\OneDrive - UMass Lowell\\FFEM_CSharp\\OutputFiles\\RVEs"); // (Filename.txt, Directory)
            ParaviewFiles paraview = new("Hex_20Pack", "C:\\Users\\Eric_Carey.UMLADCO\\OneDrive - UMass Lowell\\FFEM_CSharp\\ParaviewFiles\\Test");  // (Filename, Directory) -- DO NOT ADD .txt or .vtk to Filename
            //////////////// END INPUTS ///////////////

            // Submit model for analysis
            SubmitForAnalysis(mesh, inputfile, outputfile, paraview);
        }

        private static void SubmitForAnalysis(MeshFile mesh, InputFile inputfile, OutputFile outputfile, ParaviewFiles paraview)
        {
            // Read in mesh and Input file:
            Assembly assembly = new();
            Solver solver = new();
            mesh.ReadMeshFile(assembly);
            inputfile.ReadInputFile(ref assembly, ref solver, ref paraview);

            // Apply loading and post process:
            assembly.AnalyzeRVE(solver, outputfile, paraview);
        }
    }
}

