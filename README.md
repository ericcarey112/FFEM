# FFEM
FFEM C# source code 

## RUNNING A MODEL
So far I have a hex pack model Jamal made as an example file to run/debug the code with. Here are some steps and information to run model with what I currently have:
The main program (Which is lame) is located in repos/ConsoleApp1/Program.cs. For now, to run an example RVE, there are 4 inputs:
    1. Mesh file  (This is read in from a file. I have an example in the MeshFiles folder called Hex_20Pack.txt) <--- Currently this is generated from MATLAB
    2. Input file (This is read in from a file. I have an example Input file folder called ExampleInputFile.txt)
    3. Output file (This is generated from the program. you just need to supply a name) <--- This is where the strain/stress/stiffness for each load step is written to
    4. Paraview files (This is generated from the program. you just need to supply a name) <--- If you want to generate paraview files (Optional)
    
    Each input takes the following format: (Filename.txt, FileDirectory) with the exception being the Paraview files input which has the format (Filename, FileDirectory). There is no .txt on the end of the paraview filename.
    
    I have the file names already written, you will just need to change the file directory for each input on your computer. I'm not really sure of a better way to make these inputs? Seems a little tedious the way I have it now, but if you have better suggestions for how to organize this let me know. Thanks!

## CURRENT INPUT FILE COMMENTS
I have the current Input file (Which I will be changing to better resemble one of your JED input files with comments, better formatting, etc) set up the following way:

**FIBERS
ConsitutiveModel: 2 # 1 = Isotropic, 2 = Transversly Isotropic, 3 = Smeared Crack
E1: 276000.0 # E1 (fiber direction) in MPa
E2: 14000.0  # E2 (transverse to fiber direction) in MPa
v23: 0.25    # In-Plane poisson ratio 
G23: 20000.0 # In-Plane shear modulus in MPa

**MATRIX
ConsitutiveModel: 3
E2: 14000.0  # E2 (transverse to fiber direction) in MPa
v23: 0.25    # In-Plane poisson ratio 
Strength: 100.0 # Strength in MPa
GIC: 10.0    # GIC (Arbitrary for now)

**RVE
Thickness: 1.0 # Thickness in fiber direction
PlaneStressPlaneStrain: 2 # 1 = Plane Stress, 2 = Plane Strain

**BOUNDARY CONDITIONS
E22: 0.02 # Total applied strain in 2-direction
E33: 0.00 # Total applied strain in 3-direction
E23: 0.00 # Total applied shear strain in 2-3 plane
LoadStep: 0.00025 # Incremental strain applied for each load step

**SOLVER
Type: 2   # 1 = Static Linear, 2 = Static non-Linear
NRMax: 10 # Max number of newton-raphson iterations when searching for equilibrium
MaxAttempts: 10 # Max number of times the load step is cut in half when newton-raphson doesn't find eqbm. Simulation will end if this condition is met.
MaxLoadSteps: 500 # Max number of load steps (Either the target load E22 is reached or this MaxLoadSteps condition is reached) Simulation will end if this condition is met.
ConvergenceTolerance: 0.0025 # Force equilibrium tolerance for non-linear solutions

**PARAVIEW FILES
ParaviewVTKs: 1 # 0 = don't generate paraview vtk files, 1 = do generate paraview vtk files

## GENERAL COMMENTS:
There are still some things I know I need to change. Here is a list of immediate changes I will want to address:
    1. Change input file to include the comments I posted above and change the inputs to align more with your JED input file
    2. There are some properties and methods (Like those in abstract classes) I want to change from public to protected
    3. I have some constructors for abstract classes (Like Elements and MaterialModel) that have dummy inputs that I want to change. I have some specific questions about this, but I figured I would leave a comment about it here. 
    4. Right now all of my files are put under one project basically (FEMAssembly) because I was having problems with circular references and not really knowing how namespaces worked. I want to talk to you about how to organize these files better, but for now I left all of the "Core" files in the FFEMAssembly folder.
    5. Address more of the comments you posted on my powerpoint slides where I outlined my classes and methods
I will make changes as you post comments/suggestions whenever you have the time.

Lastly, one important thing I'm noticing is the model here in C# takes much longer to run than in MATLAB because of the matrix inversions and solving system of equations (MatrixMath.InvertMatrix and MatrixMath.LinSolve). Each load step for my example hex pack RVE takes about 5-6 seconds whereas each load step in my MATLAB model takes a split second (In the linear elastic regime). Do your load steps for FDEM run very quickly?