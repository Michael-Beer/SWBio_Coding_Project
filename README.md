# SWBio_Coding_Project
### The two example files can be found in the 'master' branch of this repository.
This is a script written in Python that analyses a Molecular Dynamics Trajectory and outputs an RMSD, RMSF, End-to-End distance and Radius of Gyration analyses into separate PNG files. This script has been run using a combination of a DCD trajectory file (NAMD, X-PLOR and CHARMM) and PDB file, but in theory should work with other trajectory files (e.g. XTC and multi-frame PDB files) and other structure files (e.g. GRO).

How To use:

1. Ensure that the computer being used to run the script has Python3.8 installed, along with MDAnalysis, numpy, matplotlib and pandas packages.

2. Place example files ('BlaC_dry.pdb' and 'BlaC.dcd') or your own structure and trajectory files into the same directory as the script.

3. If using your own structure and trajectory files, open the script and append the script to include your own file names, the residue numbers for definition of Calphas, and the residue numbers for FULL definition, and the residue number for Cterm definition - this is required so that the analyses take into account the residues of interest for your trajectory analyses. The default setting is for the BlaC.dcd trajectory.

4. Run the Python script using ./SWBIOPYTHON.py

The output of the script should be four separate PDF files named RMSD.png, RMSF.png, E2E.png and RadGyr.png containing images of the RSMD, RMSF, End-to-End and Radius of Gyration plots of your MD trajectory respectively. These will be output into the current directory. For a 1000 frame MD trajectory of a 265 amino acid protein on a decent CPU, the script should take less than 20 seconds to run.
