using myMath;

namespace FEMAssembly
{
    /// <summary>
    /// Class to generate paraview VTK files
    /// </summary>
    public class ParaviewFiles
    {
        // Properties
        public string FileName { get; set; }
        public string FileDirectory { get; set; }
        public bool GenerateFiles { get; set; }
        public int NumberOfCells { get; set; }
        public List<int> QuadElements { get; set; }
        public List<double[]> NodalLocations { get; set; }
        public List<double[]> NewDisplacements { get; set; }
        public List<int[]> Connectivity { get; set; }
        public List<double[]> NodalDisplacements { get; set; }
        public List<double[]> Strain { get; set; }
        public List<double[]> Stress { get; set; }
        public List<double[]> PrincipalStrain { get; set; }
        public List<double[]> PrincipalStress { get; set; }

        // Constructor
        public ParaviewFiles(string filename, string directory)
        {

            this.FileName = filename;
            this.FileDirectory = @directory;
            this.GenerateFiles = false;
            this.NumberOfCells = 0;
            this.QuadElements = [];
            this.NodalLocations = [];
            this.NewDisplacements = [];
            this.Connectivity = new List<int[]>();
            this.NodalDisplacements = [];
            this.Strain = [];
            this.Stress = [];
            this.PrincipalStrain = [];
            this.PrincipalStress = [];
        }

        // Methods
        /// <summary>
        /// Generate the vtk files for paraview
        /// </summary>
        public void GenerateVTKFiles(Assembly assembly, int LoadStepNumber)
        {
            // Generate Paraview files
            if (LoadStepNumber == 1) { DiscretizeCellsAroundIPs(assembly); }
            InterpolateAtNodes(assembly);
            GenerateNodalDataVTKFile(assembly, LoadStepNumber);
            DisplacementsNewCells(assembly);
            GenerateCellIPDataVTKFile(assembly, LoadStepNumber);
        }

        /// <summary>
        /// VTK files for nodal data
        /// </summary>
        private void GenerateNodalDataVTKFile(Assembly assembly, int LoadStep)
        {
            // Inputs
            List<double[]> NodalLocations = assembly.NodalLocations;
            List<int[]> ConnectMat = assembly.Connectivity;

            string filename = this.FileName + "_NodalData_" + LoadStep.ToString() + ".vtk";
            string FilePath = Path.Combine(this.FileDirectory, filename);
            using (StreamWriter writer = new StreamWriter(FilePath))
            {
                // FILE HEADER
                writer.WriteLine("# vtk DataFile Version 2.0");
                writer.WriteLine("Unstructured Grid Example");
                writer.WriteLine("ASCII");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");
                writer.WriteLine();

                // POINTS
                writer.WriteLine("POINTS " + assembly.NumberOfNodes + " float");
                for (int i = 0; i < assembly.NumberOfNodes; i++)
                {
                    double x = NodalLocations[i][0];
                    double y = NodalLocations[i][1];
                    double z = 0.0;
                    writer.WriteLine(x + " " + y + " " + z);
                }
                writer.WriteLine();

                // CELL CONNECTIVITY
                int NumLocalNodes = 0;
                for (int i = 0; i < assembly.NumberOfElements; i++)
                {
                    int[] ConnectRow = assembly.Connectivity[i];
                    for (int j = 0; j < ConnectRow.Length; j++)
                    {
                        if (ConnectRow[j] == 0) { continue; }
                        else { NumLocalNodes++; }
                    }
                }
                int TotalCellData = assembly.NumberOfElements + NumLocalNodes;
                writer.WriteLine("CELLS " + assembly.NumberOfElements + " " + TotalCellData);

                for (int i = 0; i < assembly.NumberOfElements; i++)
                {
                    int[] ConnectRow = assembly.Connectivity[i];
                    double[] NewRow = new double[ConnectRow.Length];
                    if (assembly.ElementList[i].Type == "8NQ")
                    {
                        NewRow[0] = ConnectRow[0] - 1; NewRow[1] = ConnectRow[2] - 1;
                        NewRow[2] = ConnectRow[4] - 1; NewRow[3] = ConnectRow[6] - 1;
                        NewRow[4] = ConnectRow[1] - 1; NewRow[5] = ConnectRow[3] - 1;
                        NewRow[6] = ConnectRow[5] - 1; NewRow[7] = ConnectRow[7] - 1;
                    }
                    else if (assembly.ElementList[i].Type == "6NT")
                    {
                        NewRow[0] = ConnectRow[0] - 1; NewRow[1] = ConnectRow[2] - 1;
                        NewRow[2] = ConnectRow[4] - 1; NewRow[3] = ConnectRow[1] - 1;
                        NewRow[4] = ConnectRow[3] - 1; NewRow[5] = ConnectRow[5] - 1;
                    }

                    writer.Write(assembly.ElementList[i].NumNodes);
                    for (int j = 0; j < NewRow.Length; j++)
                    {
                        writer.Write(" " + NewRow[j]);
                    }
                    writer.WriteLine();
                }
                writer.WriteLine();

                // CELL TYPES (Paraview IDs)
                int CellID = 0;
                writer.WriteLine("CELL_TYPES " + assembly.NumberOfElements);
                for (int i = 0; i < assembly.NumberOfElements; i++)
                {
                    string type = assembly.ElementList[i].Type;
                    if (type == "4NQ") { CellID = 9; }
                    else if (type == "3NT") { CellID = 5; }
                    else if (type == "6NQ") { CellID = 7; }
                    else if (type == "8NQ") { CellID = 23; }
                    else if (type == "4NT") { CellID = 7; }
                    else if (type == "6NT") { CellID = 22; }
                    else { throw new Exception("Incorrect Element Type ParaviewFiles -> GenerateNodalData"); }
                    writer.WriteLine(CellID);
                }
                writer.WriteLine();

                // POINT DATA
                writer.WriteLine("POINT_DATA " + assembly.NumberOfNodes);

                // Displacements
                writer.WriteLine("SCALARS " + "U " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, assembly.NumberOfNodes, this.NodalDisplacements);
                writer.WriteLine();

                // Strains
                writer.WriteLine("SCALARS " + "E " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, assembly.NumberOfNodes, this.Strain);
                writer.WriteLine();

                // Stresses
                writer.WriteLine("SCALARS " + "S " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, assembly.NumberOfNodes, this.Stress);
                writer.WriteLine();

                // Principal Strains
                writer.WriteLine("SCALARS " + "PE " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, assembly.NumberOfNodes, this.PrincipalStrain);
                writer.WriteLine();

                // Principal Stresses
                writer.WriteLine("SCALARS " + "PS " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, assembly.NumberOfNodes, this.PrincipalStress);
            }
        }

        /// <summary>
        /// VTK files for cell data (Damage)
        /// </summary>
        /// <param name="assembly"></param>
        /// <param name="LoadStep"></param>
        /// <exception cref="Exception"></exception>
        private void GenerateCellIPDataVTKFile(Assembly assembly, int LoadStep)
        {
            string filename = this.FileName + "_CellData_" + LoadStep.ToString() + ".vtk";
            string FilePath = Path.Combine(this.FileDirectory, filename);
            using (StreamWriter writer = new StreamWriter(FilePath))
            {
                // FILE HEADER
                writer.WriteLine("# vtk DataFile Version 2.0");
                writer.WriteLine("Unstructured Grid Example");
                writer.WriteLine("ASCII");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");
                writer.WriteLine();

                // POINTS
                writer.WriteLine("POINTS " + NodalLocations.Count + " float");
                for (int i = 0; i < NodalLocations.Count; i++)
                {
                    double x = NodalLocations[i][0];
                    double y = NodalLocations[i][1];
                    double z = 0.0;
                    writer.WriteLine(x + " " + y + " " + z);
                }
                writer.WriteLine();

                // CELL CONNECTIVITY
                int NumLocalNodes = 0;
                for (int i = 0; i < NumberOfCells; i++)
                {
                    int[] ConnectRow = this.Connectivity[i];
                    for (int j = 0; j < ConnectRow.Length; j++)
                    {
                        if (ConnectRow[j] == 0) { continue; }
                        else { NumLocalNodes++; }
                    }
                }
                int TotalCellData = NumberOfCells + NumLocalNodes;
                writer.WriteLine("CELLS " + NumberOfCells + " " + TotalCellData);

                for (int i = 0; i < NumberOfCells; i++)
                {
                    int[] ConnectRow = this.Connectivity[i];
                    double[] NewRow = new double[ConnectRow.Length];
                    if (i < assembly.NumberOfElements)
                    {
                        if (assembly.ElementList[i].Type == "8NQ")
                        {
                            NewRow[0] = ConnectRow[0] - 1; NewRow[1] = ConnectRow[2] - 1;
                            NewRow[2] = ConnectRow[4] - 1; NewRow[3] = ConnectRow[6] - 1;
                            NewRow[4] = ConnectRow[1] - 1; NewRow[5] = ConnectRow[3] - 1;
                            NewRow[6] = ConnectRow[5] - 1; NewRow[7] = ConnectRow[7] - 1;
                        }
                        else if (assembly.ElementList[i].Type == "6NT")
                        {
                            NewRow[0] = ConnectRow[0] - 1; NewRow[1] = ConnectRow[2] - 1;
                            NewRow[2] = ConnectRow[4] - 1; NewRow[3] = ConnectRow[1] - 1;
                            NewRow[4] = ConnectRow[3] - 1; NewRow[5] = ConnectRow[5] - 1;
                        }

                        writer.Write(assembly.ElementList[i].NumNodes);
                        for (int j = 0; j < NewRow.Length; j++)
                        {
                            writer.Write(" " + NewRow[j]);
                        }
                        writer.WriteLine();
                    }
                    else if (i >= assembly.NumberOfNodes)
                    {
                        NewRow[0] = ConnectRow[0] - 1; NewRow[1] = ConnectRow[1] - 1;
                        NewRow[2] = ConnectRow[2] - 1; NewRow[3] = ConnectRow[3] - 1;

                        writer.Write(4);
                        for (int j = 0; j < 4; j++)
                        {
                            writer.Write(" " + NewRow[j]);
                        }
                        writer.WriteLine();
                    }
                }
                writer.WriteLine();

                // CELL TYPES (Paraview IDs)
                int CellID = 0;
                writer.WriteLine("CELL_TYPES " + NumberOfCells);
                for (int i = 0; i < NumberOfCells; i++)
                {
                    if (i < assembly.NumberOfElements)
                    {
                        string type = assembly.ElementList[i].Type;
                        if (type == "4NQ") { CellID = 9; }
                        else if (type == "3NT") { CellID = 5; }
                        else if (type == "6NQ") { CellID = 7; }
                        else if (type == "8NQ") { CellID = 23; }
                        else if (type == "4NT") { CellID = 7; }
                        else if (type == "6NT") { CellID = 22; }
                        else { throw new Exception("Incorrect Element Type ParaviewFiles -> GenerateIPData"); }
                        writer.WriteLine(CellID);
                    }
                    else if (i >= assembly.NumberOfElements)
                    {
                        CellID = 9;
                        writer.WriteLine(CellID);
                    }
                }
                writer.WriteLine();

                // POINT DATA
                writer.WriteLine("POINT_DATA " + NodalLocations.Count);

                // Displacements
                writer.WriteLine("SCALARS " + "U " + "float64 " + 3);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                WritePointData(writer, NodalLocations.Count, NodalDisplacements);
                writer.WriteLine();

                // CELL DATA
                writer.WriteLine("CELL_DATA " + NumberOfCells);
                writer.WriteLine("SCALARS " + "Damage " + "float64 " + 1);
                writer.WriteLine("LOOKUP_TABLE " + "default");
                for (int i = 0; i < assembly.NumberOfElements; i++)
                {
                    double damage = 0.0;
                    string type = assembly.ElementList[i].Type;
                    if (type == "3NT" || type == "6NT" && assembly.ElementList[i].ModelID == 3)
                    {
                        damage = assembly.ElementList[i].Model.StateVars[1, 0];
                    }
                    writer.WriteLine(damage);
                }

                for (int i = 0; i < QuadElements.Count; i++)
                {
                    int Quad = QuadElements[i];
                    int NumIPs = assembly.ElementList[Quad - 1].NumIPs;
                    for (int j = 0; j < NumIPs; j++)
                    {
                        if (assembly.ElementList[Quad - 1].ModelID == 3)
                        {
                            double damage = assembly.ElementList[Quad - 1].Model.StateVars[1, j];
                            writer.WriteLine(damage);
                        }
                        else { writer.WriteLine(0.0); }
                    }
                }
            }
        }

        /// <summary>
        /// Write point data (Displacements, stresses, strains)
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="NumNodes"></param>
        /// <param name="Data"></param>
        private static void WritePointData(StreamWriter writer, int NumNodes, List<double[]> Data)
        {
            for (int i = 0; i < NumNodes; i++)
            {
                for (int j = 0; j < Data[i].Length; j++)
                {
                    if (j == Data[i].Length - 1) { writer.Write(Data[i][j]); }
                    else { writer.Write(Data[i][j] + " "); }
                }
                writer.WriteLine();
            }
        }

        /// <summary>
        /// Create 4-Noded elements (Cells) around integration points
        /// </summary>
        /// <param name="assembly"></param>
        public void DiscretizeCellsAroundIPs(Assembly assembly)
        {
            // Inputs
            int CellCount = assembly.NumberOfElements;
            this.NodalLocations = new List<double[]>(assembly.NodalLocations);
            this.Connectivity = assembly.Connectivity;

            FindQuadElements(assembly); // store the quad element numbers

            for (int i = 0; i < QuadElements.Count; i++)
            {
                int quad = QuadElements[i]; // Quad element number
                Elements element = assembly.ElementList[quad - 1]; // Quad Element (Base 0 indexing)
                int NumIPs = element.NumIPs; // Number of IPs per quad
                List<double[]> XI = PointsAroundIPs(NumIPs); // Get (xi,eta) points around IPs
                List<List<double[]>> Cells = CellsAroundIPs(NumIPs, XI);

                foreach (var array in XI)
                {
                    double xi = array[0]; // xi
                    double eta = array[1]; // eta

                    // Calculate (x,y) coordinate:
                    double[] N = ShapeFunctions.EvaluateShapeFunctions(element.Type, xi, eta);
                    double[,] NMat = ShapeFunctions.AssembleShapeFunctions(element.Type, N);
                    double[] xy = MatrixMath.Multiply(NMat, element.NodalLocations);

                    for (int k = 0; k < NodalLocations.Count; k++)
                    {
                        // If xy does not overlap with any existing point, then add this as a new point
                        double[] point = [NodalLocations[k][0], NodalLocations[k][1]]; // point to check for overlap
                        bool Overlap = OverlapCheck(xy, point, 0.000001);
                        if (Overlap) { break; }
                        else if (k == NodalLocations.Count - 1)
                        {
                            this.NodalLocations.Add(xy);
                        }
                    }

                }

                // Add new points to connectivity
                foreach (var cell in Cells)
                {
                    CellCount++; // Update cell counter
                    int[] row = new int[4];
                    Array.Fill(row, 0);
                    this.Connectivity.Add(row);

                    for (int k = 0; (k < cell.Count); k++)
                    {
                        double xi = cell[k][0];
                        double eta = cell[k][1];
                        double[] N = ShapeFunctions.EvaluateShapeFunctions(element.Type, xi, eta);
                        double[,] NMat = ShapeFunctions.AssembleShapeFunctions(element.Type, N);
                        double[] xy = MatrixMath.Multiply(NMat, element.NodalLocations);

                        // Assign (x,y) to new connectivity
                        for (int l = 0; l < NodalLocations.Count; l++)
                        {
                            double[] point = [NodalLocations[l][0], NodalLocations[l][1]];


                            // Check for overlapping nodes and add to connectivity
                            bool Overlap = OverlapCheck(xy, point, 0.000001);
                            if (Overlap)
                            {
                                int NodeNum = l + 1;
                                this.Connectivity[CellCount - 1][k] = NodeNum;
                            }
                        }
                    }
                }
                this.NumberOfCells = CellCount;
            }
        }

        /// <summary>
        /// Interpolate Displacements, Stresses/Strains, Principal Stresses/Strains at nodes for paraview
        /// </summary>
        private void InterpolateAtNodes(Assembly assembly)
        {
            //Inputs
            this.NodalDisplacements = Doubles.ListOfZeros(assembly.NumberOfNodes, 3);
            this.Strain = Doubles.ListOfZeros(assembly.NumberOfNodes, 3);
            this.Stress = Doubles.ListOfZeros(assembly.NumberOfNodes, 3);
            this.PrincipalStrain = Doubles.ListOfZeros(assembly.NumberOfNodes, 3);
            this.PrincipalStress = Doubles.ListOfZeros(assembly.NumberOfNodes, 3);
            bool RVE = false;

            // Loop over elements and interpolate data at nodes
            for (int i = 0; i < assembly.NumberOfElements; i++)
            {
                Elements element = assembly.ElementList[i];
                string type = element.Type;
                int NumNodes = element.NumNodes;
                double[] NodalLocations = element.NodalLocations;
                List<int[]> ConnectMat = assembly.Connectivity;

                // Distribute displacements to nodes:
                double[] q = Assembly.GlobalToLocalVector(element, assembly.Connectivity, i, assembly.GlobalQ);

                // If the model is an RVE and the element type is a quad
                if (RVE)
                {
                    if (type == "6NQ" || type == "8NQ")
                    {
                        InterpolateQuadsForRVE(assembly, element);
                        return;
                    }
                }

                List<double[]> XiEta = GetXiEtaNodes(type, NumNodes);
                for (int j = 0; j < NumNodes; j++)
                {
                    int NodeNum = ConnectMat[i][j]; // Global node number
                    double xi = XiEta[j][0];
                    double eta = XiEta[j][1];

                    // Calculate Displacements:
                    double[] N = ShapeFunctions.EvaluateShapeFunctions(element.Type, xi, eta);
                    double[,] NMat = ShapeFunctions.AssembleShapeFunctions(element.Type, N);
                    double[] disp = MatrixMath.Multiply(NMat, q);

                    // Calculate B Matrix:
                    double[,] BMatrix = Elements.BMatrix2D(element.Type, element.NodalLocations, xi, eta);

                    // Calculate Stress/Strain:
                    double[,] DMatrix = new double[3, 3];
                    // Isotropic
                    if (element.ModelID == 1)
                    {
                        DMatrix = IsotropicLinearElastic.CalcDMatrixIsotropic(element.Model.E2, element.Model.nu23, element.PlaneStressPlaneStrain);
                    }
                    // Transversly isotropic
                    else if (element.ModelID == 2)
                    {
                        DMatrix = TransverslyIsotropicLinearElastic.CalcDMatrixTransverslyIsotropic(element.Model.E2, element.Model.nu23, element.Model.G23);
                    }
                    // Smeared Crack
                    else if (element.ModelID == 3)
                    {
                        // NEED TO CHANGE HOW IM FINDING DAMAGE. CURRENTLY IT IS SET TO 0.0 UNTIL I FIND A WAY TO GET THE RIGHT DAMAGE
                        DMatrix = SmearedCrack.CalcDMatrixIsotropicDamage(element.Model.E2, element.Model.nu23, element.PlaneStressPlaneStrain, element.Model.StateVars[1, 0]);
                    }

                    double[] Eps = Elements.CalcStrain(BMatrix, q);
                    double[] Sig = Elements.CalcStress(DMatrix, Eps);

                    // Calculate Principal Stress/Strain at Nodes:
                    double NormToCrackAngle = CalcNormToCrackAngle_Nodes(Sig);
                    double[] Eps_P = Elements.RotateStressOrStrain(Eps, NormToCrackAngle);
                    double[] Sig_P = Elements.RotateStressOrStrain(Sig, NormToCrackAngle);

                    // Store data:
                    this.NodalDisplacements[NodeNum - 1] = [disp[0], disp[1], 0.0];
                    this.Strain[NodeNum - 1] = [Eps[0], Eps[1], Eps[2]];
                    this.Stress[NodeNum - 1] = [Sig[0], Sig[1], Sig[2]];
                    this.PrincipalStrain[NodeNum - 1] = [Eps_P[0], Eps_P[1], Eps_P[2]];
                    this.PrincipalStress[NodeNum - 1] = [Sig_P[0], Sig_P[1], Sig_P[2]];
                }
            }
        }

        /// <summary>
        /// Calcualates the displacements for the new cells created for paraview
        /// </summary>
        private void DisplacementsNewCells(Assembly assembly)
        {
            // Calculate displacements around for new cells
            int NumElemsPlusQuads = assembly.NumberOfElements + QuadElements.Count;
            this.NodalDisplacements = Doubles.ListOfZeros(this.NodalLocations.Count, 3);
            int QuadCount = 0;
            int RowCount = 0;
            for (int i = 0; i < NumElemsPlusQuads; i++)
            {
                if (i < assembly.NumberOfElements - 1)
                {
                    int NumNodes = assembly.ElementList[i].NumNodes;
                    for (int j = 0; j < NumNodes; j++)
                    {
                        int GlobalNode = this.Connectivity[i][j];
                        int xDOF = assembly.NDOFPNode * GlobalNode - 2;
                        int yDOF = assembly.NDOFPNode * GlobalNode - 1;
                        this.NodalDisplacements[GlobalNode - 1][0] = assembly.GlobalQ[xDOF];
                        this.NodalDisplacements[GlobalNode - 1][1] = assembly.GlobalQ[yDOF];
                        this.NodalDisplacements[GlobalNode - 1][2] = 0.0;
                    }
                }
                else if (i >= assembly.NumberOfElements)
                {
                    // Get local dispalcements for new cell
                    QuadCount++;
                    int quad = QuadElements[QuadCount - 1];
                    Elements element = assembly.ElementList[quad - 1];
                    int NumIPs = element.NumIPs;
                    double[] q = Assembly.GlobalToLocalVector(element, assembly.Connectivity, quad - 1, assembly.GlobalQ);

                    // Get (xi, eta) points for new cell
                    List<double[]> XI = PointsAroundIPs(NumIPs);
                    List<List<double[]>> Cells = CellsAroundIPs(NumIPs, XI);

                    int CellsInElement = Cells.Count;
                    for (int j = 0; j < CellsInElement; j++)
                    {
                        RowCount++;
                        int RowNum = assembly.NumberOfElements + RowCount;
                        for (int k = 0; k < Cells[j].Count; k++)
                        {
                            int NodeNum = this.Connectivity[RowNum - 1][k];

                            // (xi, eta) at corner of cell
                            double xi = Cells[j][k][0];
                            double eta = Cells[j][k][1];

                            // Calculate displacements
                            double[] N = ShapeFunctions.EvaluateShapeFunctions(element.Type, xi, eta);
                            double[,] NMat = ShapeFunctions.AssembleShapeFunctions(element.Type, N);
                            double[] disp = MatrixMath.Multiply(NMat, q);
                            this.NodalDisplacements[NodeNum - 1][0] = disp[0];
                            this.NodalDisplacements[NodeNum - 1][1] = disp[1];
                            this.NodalDisplacements[NodeNum - 1][2] = 0.0;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Interpolates at nodes specifically for 8-Noded quads (Most nodes are shared with fiber elements)
        /// </summary>
        private void InterpolateQuadsForRVE(Assembly assembly, Elements element)
        {
            List<double[]> XiEta = GetXiEtaNodes(element.Type, element.NumNodes);
            for (int i = 0; i < element.NumNodes; i++)
            {
                int NodeNum = (int)assembly.Connectivity[element.ElementNumber - 1][i]; // Global node number

                // Check if data has already been calculated at this node
                double DataCheck = this.NodalDisplacements[NodeNum - 1][0];
                if (DataCheck == 0.0)
                {
                    double xi = XiEta[i][0];
                    double eta = XiEta[i][1];

                    // Calculate Displacements:
                    double[] N = ShapeFunctions.EvaluateShapeFunctions(element.Type, xi, eta);
                    double[,] NMat = ShapeFunctions.AssembleShapeFunctions(element.Type, N);
                    double[] disp = MatrixMath.Multiply(NMat, element.NodalLocations);

                    // Determine if nodes are adjecent to midpoint nodes:
                    int NodeNum_L;
                    int NodeNum_R;
                    if (i == 7)
                    {
                        NodeNum_L = assembly.Connectivity[element.ElementNumber - 1][6];
                        NodeNum_R = assembly.Connectivity[element.ElementNumber - 1][0];
                    }
                    else if (i == 0)
                    {
                        NodeNum_L = assembly.Connectivity[element.ElementNumber - 1][7];
                        NodeNum_R = assembly.Connectivity[element.ElementNumber - 1][1];
                    }
                    else
                    {
                        NodeNum_L = (int)assembly.Connectivity[element.ElementNumber - 1][i - 1];
                        NodeNum_R = (int)assembly.Connectivity[element.ElementNumber - 1][i + 1];
                    }

                    // Calculate strains/stresses and principal strains/stresses:
                    double[] Eps = new double[3];
                    double[] Sig = new double[3];
                    double[] Eps_P = new double[3];
                    double[] Sig_P = new double[3];
                    Eps[0] = (this.Strain[NodeNum_L - 1][0] + this.Strain[NodeNum_R - 1][0]) / 2.0;
                    Eps[1] = (this.Strain[NodeNum_L - 1][1] + this.Strain[NodeNum_R - 1][1]) / 2.0;
                    Eps[2] = (this.Strain[NodeNum_L - 1][2] + this.Strain[NodeNum_R - 1][2]) / 2.0;
                    Sig[0] = (this.Stress[NodeNum_L - 1][0] + this.Stress[NodeNum_R - 1][0]) / 2.0;
                    Sig[0] = (this.Stress[NodeNum_L - 1][1] + this.Stress[NodeNum_R - 1][1]) / 2.0;
                    Sig[0] = (this.Stress[NodeNum_L - 1][2] + this.Stress[NodeNum_R - 1][2]) / 2.0;
                    Eps_P[0] = (this.PrincipalStrain[NodeNum_L - 1][0] + this.PrincipalStrain[NodeNum_R - 1][0]) / 2.0;
                    Eps_P[1] = (this.PrincipalStrain[NodeNum_L - 1][1] + this.PrincipalStrain[NodeNum_R - 1][1]) / 2.0;
                    Eps_P[2] = (this.PrincipalStrain[NodeNum_L - 1][2] + this.PrincipalStrain[NodeNum_R - 1][2]) / 2.0;
                    Sig_P[0] = (this.PrincipalStress[NodeNum_L - 1][0] + this.PrincipalStress[NodeNum_R - 1][0]) / 2.0;
                    Sig_P[1] = (this.PrincipalStress[NodeNum_L - 1][1] + this.PrincipalStress[NodeNum_R - 1][1]) / 2.0;
                    Sig_P[2] = (this.PrincipalStress[NodeNum_L - 1][2] + this.PrincipalStress[NodeNum_R - 1][2]) / 2.0;

                    // Store data:
                    this.NodalDisplacements[NodeNum - 1] = [disp[0], disp[1], 0.0];
                    this.Strain[NodeNum - 1] = [Eps[0], Eps[1], Eps[2]];
                    this.Stress[NodeNum - 1] = [Sig[0], Sig[1], Sig[2]];
                    this.PrincipalStrain[NodeNum - 1] = [Eps_P[0], Eps_P[1], Eps_P[2]];
                    this.PrincipalStress[NodeNum - 1] = [Sig_P[0], Sig_P[1], Sig_P[2]];
                }
            }
        }

        /// <summary>
        /// Obtain (xi,eta) coordinates of nodes based on element type and number of nodes
        /// </summary>
        /// <param name="type"></param>
        /// <param name="NumNodes"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private static List<double[]> GetXiEtaNodes(string type, int NumNodes)
        {
            List<double[]> XiEta = [];
            if (type == "3NT" || type == "4NT" || type == "6NT")
            {
                switch (NumNodes)
                {
                    case 3:
                        {
                            XiEta.Add([0, 0]); XiEta.Add([1, 0]); XiEta.Add([0, 1]);
                            break;
                        }
                    case 4:
                        {
                            XiEta.Add([0, 0.99999]); XiEta.Add([0, 0]);
                            XiEta.Add([0.5, 0]); XiEta.Add([1, 0]);
                            break;
                        }
                    case 6:
                        {
                            XiEta.Add([0, 1]); XiEta.Add([0, 0.5]); XiEta.Add([0, 0]);
                            XiEta.Add([0.5, 0]); XiEta.Add([1, 0]); XiEta.Add([0.5, 0.5]);
                            break;
                        }
                    default:
                        {
                            throw new Exception("Incorrect number of IPs for triangle element Paraview -> XiEtaNodes");
                        }
                }
            }
            else if (type == "4NQ" || type == "6NQ" || type == "8NQ")
            {
                switch (NumNodes)
                {
                    case 4:
                        {
                            XiEta.Add([-1, -1]); XiEta.Add([1, -1]);
                            XiEta.Add([1, 1]); XiEta.Add([-1, 1]);
                            break;
                        }
                    case 6:
                        {
                            XiEta.Add([-1, -1]); XiEta.Add([0, -1]); XiEta.Add([1, -1]);
                            XiEta.Add([1, 1]); XiEta.Add([0, 1]); XiEta.Add([-1, 1]);
                            break;
                        }
                    case 8:
                        {
                            XiEta.Add([-1, -1]); XiEta.Add([0, -1]); XiEta.Add([1, -1]); XiEta.Add([1, 0]);
                            XiEta.Add([1, 1]); XiEta.Add([0, 1]); XiEta.Add([-1, 1]); XiEta.Add([-1, 0]);
                            break;
                        }
                    default:
                        {
                            throw new Exception("Incorrect number of IPs for quad element Paraview -> XiEtaNodes");
                        }
                }
            }
            return XiEta;
        }

        /// <summary>
        /// Finds all quad element numbers
        /// </summary>
        private void FindQuadElements(Assembly assembly)
        {
            // Loop through all elements and find quads:
            for (int i = 0; i < assembly.NumberOfElements; i++)
            {
                // Get element type:
                string type = assembly.ElementList[i].Type;
                if (type == "4NQ" || type == "6NQ" || type == "8NQ")
                {
                    this.QuadElements.Add(i + 1); // Store element number
                }
                else { continue; }
            }
        }

        /// <summary>
        /// (xi,eta) coordinates for cells around each IP
        /// </summary>
        /// <param name="NumIPs"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private static List<double[]> PointsAroundIPs(int NumIPs)
        {
            List<double[]> XI = new List<double[]>();
            switch (NumIPs)
            {
                case 1:
                    {
                        double[] array1 = [-1, -1]; double[] array2 = [1, -1];
                        double[] array3 = [1, 1]; double[] array4 = [-1, 1];
                        XI.Add(array1); XI.Add(array2);
                        XI.Add(array3); XI.Add(array4);
                        break;
                    }
                case 4:
                    {
                        double[] array1 = [-1, -1]; double[] array2 = [0, -1]; double[] array3 = [1, -1];
                        double[] array4 = [-1, 0]; double[] array5 = [0, 0]; double[] array6 = [1, 0];
                        double[] array7 = [-1, 1]; double[] array8 = [0, 1]; double[] array9 = [1, 1];
                        XI.Add(array1); XI.Add(array2); XI.Add(array3);
                        XI.Add(array4); XI.Add(array5); XI.Add(array6);
                        XI.Add(array7); XI.Add(array8); XI.Add(array9);
                        break;
                    }
                case 9:
                    {
                        double[] array1 = [-1, -1]; double[] array2 = [-1.0 / 3.0, -1]; double[] array3 = [1.0 / 3.0, -1]; double[] array4 = [1, -1];
                        double[] array5 = [-1, -1.0 / 3.0]; double[] array6 = [-1.0 / 3.0, -1.0 / 3.0]; double[] array7 = [1.0 / 3.0, -1.0 / 3.0]; double[] array8 = [1, -1.0 / 3.0];
                        double[] array9 = [-1, 1.0 / 3.0]; double[] array10 = [-1.0 / 3.0, 1.0 / 3.0]; double[] array11 = [1.0 / 3.0, 1.0 / 3.0]; double[] array12 = [1, 1.0 / 3.0];
                        double[] array13 = [-1, 1]; double[] array14 = [-1.0 / 3.0, 1]; double[] array15 = [1.0 / 3.0, 1]; double[] array16 = [1.0, 1.0];
                        XI.Add(array1); XI.Add(array2); XI.Add(array3); XI.Add(array4);
                        XI.Add(array5); XI.Add(array6); XI.Add(array7); XI.Add(array8);
                        XI.Add(array9); XI.Add(array10); XI.Add(array11); XI.Add(array12);
                        XI.Add(array13); XI.Add(array14); XI.Add(array15); XI.Add(array16);
                        break;
                    }
                default:
                    {
                        throw new Exception("Incorrect number of IPs inside Paraview>PointsAroundIPs");
                    }

            }
            return XI;
        }

        /// <summary>
        /// Add the (xi,eta) coordinates to create cells around the IPs
        /// </summary>
        /// <param name="NumIPs"></param>
        /// <param name="XI"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        private static List<List<double[]>> CellsAroundIPs(int NumIPs, List<double[]> XI)
        {
            List<List<double[]>> Cells = [];

            switch (NumIPs)
            {
                case 1:
                    {
                        Cells.Add(new List<double[]> { XI[0], XI[1], XI[2], XI[3] });
                        break;
                    }
                case 4:
                    {
                        Cells.Add(new List<double[]> { XI[4], XI[5], XI[8], XI[7] });
                        Cells.Add(new List<double[]> { XI[1], XI[2], XI[5], XI[4] });
                        Cells.Add(new List<double[]> { XI[3], XI[4], XI[7], XI[6] });
                        Cells.Add(new List<double[]> { XI[0], XI[1], XI[4], XI[3] });
                        break;
                    }
                case 9:
                    {
                        Cells.Add(new List<double[]> { XI[10], XI[11], XI[15], XI[14] });
                        Cells.Add(new List<double[]> { XI[6], XI[7], XI[11], XI[10] });
                        Cells.Add(new List<double[]> { XI[2], XI[3], XI[7], XI[6] });
                        Cells.Add(new List<double[]> { XI[9], XI[10], XI[14], XI[13] });
                        Cells.Add(new List<double[]> { XI[5], XI[6], XI[10], XI[9] });
                        Cells.Add(new List<double[]> { XI[1], XI[2], XI[6], XI[5] });
                        Cells.Add(new List<double[]> { XI[8], XI[9], XI[13], XI[12] });
                        Cells.Add(new List<double[]> { XI[4], XI[5], XI[9], XI[8] });
                        Cells.Add(new List<double[]> { XI[0], XI[1], XI[5], XI[4] });
                        break;
                    }
                default:
                    {
                        throw new Exception("Incorrect number of IPs inside Paraview>PointsAroundIPs");
                    }
            }
            return Cells;
        }

        /// <summary>
        /// Checks for overlap between 2 points
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="delta"></param>
        /// <returns></returns>
        private static bool OverlapCheck(double[] p1, double[] p2, double delta)
        {
            // Inputs
            double x1 = p1[0]; double x2 = p2[0];
            double y1 = p1[1]; double y2 = p2[1];

            double dx = Math.Abs(x2 - x1);
            double dy = Math.Abs(y2 - y1);

            // Set Overlap to true if there is overlap
            bool Overlap = false;
            if (dx <= delta && dy <= delta)
            {
                Overlap = true;
            }

            return Overlap;
        }

        /// <summary>
        /// Calcualates normal to crack angle at nodes (Note this is different than the NormToCrackAngle function in elements since those specifically return angles for an IP. This is for nodes for paraview)
        /// </summary>
        /// <returns></returns>
        private static double CalcNormToCrackAngle_Nodes(double[] Stress)
        {
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
    }
}
