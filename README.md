# Holoware-TestCases

_Input and structure files for the paper "Comparing the performance of different AMBER/GAFF forcefields, partial charge assignments, and water models in absolute binding free energy calculations"_

Input files are organized as follows:

- Forcefield parameters are in the "parameters" directory
- Input files are in the directories named by test case: "BRD4", "cMET", and "PDE2A". 
- Within each test case directory are the 16 parameter set directories. 
- Within each parameter set directory are the ligand directories.
- Within each ligand directory are the python input files, an input directory, an output directory, and a link to the parameters directory

To run the calculations for a given ligand, use the following commands and redirect the output to the output directory as indicated: 

- python solvent.py > output/solvent.out
- python restraint.py > output/restraint.out
- python complex.py > output/complex.out

Once all calculations are complete, run the following command:

- python results.py

The results file is named ddG.res

The python code requires openmm and openmmtools. It is known to work with the following python packages:

	openmm                    7.4.2           py36_cuda101_rc_1    http://conda.binstar.org/omnia
	openmmtools               0.19.0                   py36_0    http://conda.binstar.org/omnia

Protein and ligand structure files are in "Structures" directory.

	Protein structure : protein.pdb
	Ligand structures (AM1-BCC charges): ligands_am1bcc.mol2
	Ligand structures (RESP charges): ligands_jaguar.mol2
