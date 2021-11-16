# Holoware-TestCases

Input and structure files for the paper "Comparing the performance of different AMBER/GAFF forcefields, partial charge assignments, and water models in absolute binding free energy calculations"

Input files are in the directories named by test case. Within each test case directory are the 16 parameter set directories. Within each parameter set directory are the ligand directories.

Each ligand directory has the following files

The python code works with the following python packages:

	openmm                    7.4.2           py36_cuda101_rc_1    http://conda.binstar.org/omnia
	openmmtools               0.19.0                   py36_0    http://conda.binstar.org/omnia

Protein and ligand structure files are in "Structures" directory.

	Protein structure : protein.pdb
	Ligand structures (AM1-BCC charges): ligands_am1bcc.mol2
	Ligand structures (RESP charges): ligands_jaguar.mol2

Contact: David Huggins <dhuggins@tritdi.org>
