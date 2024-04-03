namespace FEMAssembly
{
    /// <summary>
    /// Class to generate output files
    /// </summary>
    public class OutputFile
    {
        // Properties
        public string FileName { get; set; }
        public string FileDirectory { get; set; }
        public List<double> Data { get; set; }

        // Constructor
        public OutputFile(string filename, string directory)
        {
            this.FileName = filename;
            this.FileDirectory = @directory;
            this.Data = new List<double>(36);
        }

        // Methods
        /// <summary>
        /// Creates an output file containing solver information, homogenized strain, stress, and stiffness
        /// </summary>
        public void UpdateOutputFile(Assembly assembly, Solver solver)
        {
            string FilePath = Path.Combine(this.FileDirectory, this.FileName);

            // Initialize header on first load step
            if (solver.LoadStepNumber == 1)
            {
                using (StreamWriter writer = new StreamWriter(FilePath))
                {
                    writer.Write("Step|Attempts|NR_Iter|      E11     |     E22     |     E33     |     E12     |     E13     |     E23     |     S11     |     S22     |     S33     |     S12     |     S13     |     S23     |     C11     |     C12     |     C13     |     C14     |     C15     |     C16     |     C22     |     C23     |     C24     |     C25     |     C26     |     C33     |     C34     |     C35     |     C36     |     C44     |     C45     |     C46     |     C55     |     C56     |     C66     |");
                }
            }

            // Solver info
            Data.Add(solver.LoadStepNumber);
            Data.Add(solver.AttemptCounter);
            Data.Add(solver.NRCounter);

            // Strains
            Data.Add(0.0); // E11
            Data.Add(assembly.HomogenizedStrain[0]); // E22
            Data.Add(assembly.HomogenizedStrain[1]); // E33
            Data.Add(0.0); // E12
            Data.Add(0.0); // E13
            Data.Add(assembly.HomogenizedStrain[2]); // E23

            // Stresses
            Data.Add(0.0); // S11
            Data.Add(assembly.HomogenizedStress[0]); // S22
            Data.Add(assembly.HomogenizedStress[1]); // S33
            Data.Add(0.0); // S12 
            Data.Add(0.0); // S13
            Data.Add(assembly.HomogenizedStress[2]); // S23

            // Stiffness
            Data.Add(0.0); Data.Add(0.0); Data.Add(0.0); // C11, C12, C13
            Data.Add(0.0); Data.Add(0.0); Data.Add(0.0); // C14, C15, C16
            Data.Add(assembly.HomogenizedStiffness[0, 0]); Data.Add(assembly.HomogenizedStiffness[0, 1]); Data.Add(0.0); // C22, C23, C24
            Data.Add(0.0); Data.Add(assembly.HomogenizedStiffness[0, 2]); Data.Add(assembly.HomogenizedStiffness[1, 1]); // C25, C26, C33
            Data.Add(0.0); Data.Add(0.0); Data.Add(assembly.HomogenizedStiffness[1, 2]); // C34, C35, C36
            Data.Add(0.0); Data.Add(0.0); Data.Add(0.0);   // C44, C45, C46
            Data.Add(0.0); Data.Add(0.0); Data.Add(assembly.HomogenizedStiffness[2, 2]); // C55 C56 C66


            using (StreamWriter writer = new StreamWriter(FilePath, true))
            {
                writer.WriteLine();
            }

            // Spacing is made manually so output file makes a table. Not sure best way to do this?
            for (int i = 0; i < Data.Count; i++)
            {
                string str;
                double data = Data[i];

                // Load Steps
                if (i == 0)
                {
                    if (data < 10) { str = "       "; }
                    else if (data >= 10 && data < 100) { str = "      "; }
                    else if (data >= 100 && data < 1000) { str = "     "; }
                    else if (data >= 1000 && data < 10000) { str = "    "; }
                    else { str = "   "; }

                    using (StreamWriter writer = new StreamWriter(FilePath, true))
                    {
                        writer.Write(data.ToString("0") + str);
                    }

                }

                // Attempts
                if (i == 1)
                {
                    if (data < 10) { str = "        "; }
                    else if (data >= 10 && data < 100) { str = "       "; }
                    else if (data >= 100 && data < 1000) { str = "      "; }
                    else if (data >= 1000 && data < 10000) { str = "     "; }
                    else { str = "    "; }

                    using (StreamWriter writer = new StreamWriter(FilePath, true))
                    {
                        writer.Write(data.ToString("0") + str);
                    }

                }

                // NR Iterations
                if (i == 2)
                {
                    if (data < 10) { str = "      "; }
                    else if (data >= 10 && data < 100) { str = "     "; }
                    else if (data >= 100 && data < 1000) { str = "    "; }
                    else if (data >= 1000 && data < 10000) { str = "   "; }
                    else { str = "  "; }

                    using (StreamWriter writer = new StreamWriter(FilePath, true))
                    {
                        writer.Write(data.ToString("0") + str);
                    }

                }

                // Strain/Stress/Stiffness
                if (i > 2)
                {
                    if (data == 0.0) { str = "    "; }
                    else if (data < 1.0 && data > 0.0) { str = "   "; }
                    else if (data < 0.0) { str = "   "; }
                    else { str = "    "; }

                    using (StreamWriter writer = new StreamWriter(FilePath, true))
                    {
                        writer.Write(data.ToString("0.00000E00") + str);
                    }
                }
            }
            Data.Clear();
        }
    }
}
