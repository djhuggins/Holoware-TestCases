#!/usr/bin/env python
# coding: utf-8
import math
import openmmtools
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
from sys import stdout
import numpy as np
from importlib import reload
from simtk.openmm import Vec3
from simtk import unit
from openmmtools.multistate import MultiStateReporter, MultiStateSampler, ReplicaExchangeSampler, ParallelTemperingSampler, SAMSSampler, ReplicaExchangeAnalyzer
from openmmtools.states import GlobalParameterState
import os

kB = 0.008314472471220214 * unit.kilojoules_per_mole/unit.kelvin
temp = 300*unit.kelvin
kT = kB * temp
kcal = 4.1868 * unit.kilojoules_per_mole
kTtokcal = kT/kcal * unit.kilocalories_per_mole
runflag='run'

class RestraintComposableState(GlobalParameterState):
	lambda_restraints=GlobalParameterState.GlobalParameter('lambda_restraints', standard_value=1.0)

args = sys.argv[1:]
if len(args) > 1:
	raise ValueError('Only take one argument runflag')

if len(args) == 0:
	runflag='run'
else:
	runflag=args[0]
	allowedrunflags = ['run', 'extend', 'recover']
	if runflag not in allowedrunflags:
		raise ValueError('Please select runflag from {}'.format(allowedrunflags))

#Forcefield
proteinXML='parameters/protein.fb15.xml'
forcefieldXML='parameters/gaff-2.11.xml'
waterModel='tip4pfb'

#System
ionicStrength=0.15*unit.molar
boxPadding=0.9*unit.nanometers

#Pathway
stericsSteps=24
elecSteps=12
restraintSteps=8
stericsExponent=1.3
elecExponent=0.5
restraintExponent=2.5
stericsInversion=0.15
elecInversion=0.5
restraintInversion=0

#Dynamics
timeStep=4.0*unit.femtoseconds
stepsPerIteration=500
productionIterations=1000
equilibrationIterations=200
iterationsPerCheckpoint=100
extendIterations=1000

waterXML='parameters/'+waterModel+'_standard.xml'
ionXML='parameters/'+waterModel+'_HFE_multivalent.xml'

print("Padding:", boxPadding)
print("Forcefield:", proteinXML, " ", forcefieldXML, " ", waterXML, " ", ionXML)

pdb = app.PDBFile('input/solvent.pdb')
forcefield = app.ForceField(proteinXML, waterXML, ionXML, forcefieldXML, 'input/ligand.xml')
solvated = app.Modeller (pdb.topology, pdb.positions)
solvated.addExtraParticles(forcefield)
solvated.addSolvent(forcefield, model='tip4pew', ionicStrength=ionicStrength, neutralize=True, padding=boxPadding)
system = forcefield.createSystem(solvated.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005, hydrogenMass=4*unit.amus)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))
forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
nforces = system.getNumForces()

restraint_state = RestraintComposableState(lambda_restraints=1.0)
ligand_a_atoms = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]
ligand_a_bonds = []
ligand_a_angles =[]
ligand_a_torsions =[]
bond_energy_function = "lambda_restraints;"
harmonicforce=mm.CustomBondForce(bond_energy_function)
harmonicforce.addGlobalParameter('lambda_restraints', 1.0)
system.addForce(harmonicforce)
#Get date from forcefield
bond_force=forces['HarmonicBondForce']
angle_force = forces['HarmonicAngleForce']
torsion_force=forces['PeriodicTorsionForce']
bond_data = forcefield.getGenerators()[0]
angle_data = forcefield.getGenerators()[1]
torsion_data = forcefield.getGenerators()[2]
nonbonded_data = forcefield.getGenerators()[3]
num_nonbonded=len(nonbonded_data.params.paramsForType)
num_bonds=len(bond_data.length)
num_angles=len(angle_data.angle)
num_propers=len(torsion_data.proper)
num_impropers=len(torsion_data.improper)

num_a_atoms=len(ligand_a_atoms)
num_a_bonds=len(ligand_a_bonds)
num_a_angles=len(ligand_a_angles)
num_a_torsions=len(ligand_a_torsions)

for index in range(bond_force.getNumBonds()):
        [atom_i, atom_j, r0, K] = bond_force.getBondParameters(index)
        if set([atom_i, atom_j]).intersection(ligand_a_atoms):
                if set([atom_i]).intersection(ligand_a_atoms) and not set([atom_j]).intersection(ligand_a_atoms):
                        print("Fatal Error: bond between ligand",atom_i,"and non ligand",atom_j,". Exiting.")
                        exit(0)
                if set([atom_j]).intersection(ligand_a_atoms) and not set([atom_i]).intersection(ligand_a_atoms):
                        print("Fatal Error: bond between ligand",atom_j,"and non ligand",atom_i,". Exiting.")
                        exit(0)

#Setup alchemical system
reload(openmmtools.alchemy)
factory = openmmtools.alchemy.AbsoluteAlchemicalFactory(consistent_exceptions=False, split_alchemical_forces = True, alchemical_pme_treatment = 'exact')
reference_system = system

#OpenMM crashes if one adds nonexistent alchemical angles or torsions
if(num_a_bonds == 0):
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms,alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')
else:
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')

alchemical_system_in = factory.create_alchemical_system(reference_system, alchemical_regions = [alchemical_region_A])

alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=300*unit.kelvin, pressure=1*unit.bar)
composable_states = [alchemical_state_A, restraint_state]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)
integrator=mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
context = compound_state.create_context(integrator)
alchemical_system_in=context.getSystem()
alchemical_forces = {alchemical_system_in.getForce(index).__class__.__name__: alchemical_system_in.getForce(index) for index in range(alchemical_system_in.getNumForces())}
nonbonded_force = alchemical_forces['NonbondedForce']

#Use offsets to interpolate
nonbonded_data = forcefield.getGenerators()[3]
alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=300*unit.kelvin, pressure=1*unit.bar)
composable_states = [alchemical_state_A, restraint_state]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)

#DEBUG info
sys = compound_state.get_system()
#file = open('DEBUG_solvent.xml','w')
#file.write(mm.XmlSerializer.serialize(sys))
#file.close()

stericsSteps+=1
elecSteps+=1
nstates=stericsSteps+elecSteps-1
print("There will be ", nstates, " states in total")
print("Lambda exponent: ", stericsExponent)
print("stepsPerIteration:", stepsPerIteration, " productionIterations: ", productionIterations, "equilibrationIterations: ", equilibrationIterations)
print("Timestep: ", timeStep)
box_vec = alchemical_system_in.getDefaultPeriodicBoxVectors()
print("Box vectors:", box_vec)

#Setup lambda schedule
lambdas = []
for j in range(nstates):
		column = []
		for i in range(2):
				column.append(0)
		lambdas.append(column)

lambdas_elec_A = np.linspace(0.0, 1.0, elecSteps)
lambdas_sterics_A = np.linspace(0.0, 1.0, stericsSteps)

for j in range(stericsSteps):
        if(lambdas_sterics_A[j] < stericsInversion):
                lambdas_sterics_A[j]=stericsInversion*pow(lambdas_sterics_A[j]/stericsInversion,1/stericsExponent)
        elif(lambdas_sterics_A[j] > stericsInversion):
                lambdas_sterics_A[j]=stericsInversion+(1-stericsInversion)*pow((lambdas_sterics_A[j]-stericsInversion)/(1-stericsInversion),stericsExponent)

for j in range(elecSteps):
        if(lambdas_elec_A[j] < elecInversion):
                lambdas_elec_A[j]=elecInversion*pow(lambdas_elec_A[j]/elecInversion,1/elecExponent)
        elif(lambdas_elec_A[j] > elecInversion):
                lambdas_elec_A[j]=elecInversion+(1-elecInversion)*pow((lambdas_elec_A[j]-elecInversion)/(1-elecInversion),elecExponent)

#First switch on sterics
for j in range(stericsSteps):
        lambdas[j][0]=lambdas_sterics_A[j]

#Sterics stay on
for j in range(elecSteps):
        lambdas[j+stericsSteps-1][0]=1.0

#Then elec
for j in range(elecSteps):
        lambdas[j+stericsSteps-1][1]=lambdas_elec_A[j]

#Sanity check
print("")
print("Lambdas matrix")
print("Lsterics_A      Lelec_A")
for j in range(len(lambdas)):
	for i in range(len(lambdas[j])):
		print("%-15.2f" % lambdas[j][i], end=' ')
	print("")
print("")

sampler_states = list()
thermodynamic_states = list()

for k in range(nstates):
    compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
    if(num_a_atoms != 0):
        compound_state.lambda_sterics_A=lambdas[k][0]
        compound_state.lambda_electrostatics_A=lambdas[k][1]
    if(num_a_bonds != 0):
        compound_state.lambda_bonds_A=lambdas[k][1]
    if(num_a_angles != 0):
        compound_state.lambda_angles_A=lambdas[k][1]
    if(num_a_torsions != 0):
        compound_state.lambda_torsions_A=lambdas[k][1]
    compound_state.lambda_restraints=1.0
    sys = compound_state.get_system()
    sampler_states.append(openmmtools.states.SamplerState(positions=solvated.positions, box_vectors=box_vec))
    thermodynamic_states.append(compound_state)

print("Integrator: LangevinSplittingDynamicsMove")
print("Sampler: ReplicaExchangeSampler")

lsd_move = openmmtools.mcmc.LangevinSplittingDynamicsMove(timestep=timeStep, collision_rate=1.0/unit.picoseconds, n_steps=stepsPerIteration)
print('Minimizing......')
for k in range(nstates):
	sampler_state = sampler_states[k]
	thermodynamic_state = thermodynamic_states[k]
	integrator=mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
	context = thermodynamic_state.create_context(integrator)
	system = context.getSystem()
	for force in system.getForces():
		if isinstance(force, mm.CustomBondForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.HarmonicBondForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.HarmonicAngleForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.PeriodicTorsionForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomAngleForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.NonbondedForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomNonbondedForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomTorsionForce):
			force.updateParametersInContext(context)
	sampler_state.apply_to_context(context)
	initial_energy = thermodynamic_state.reduced_potential(context)
	print("Sampler state {}: initial energy {:8.3f}kT".format(k, initial_energy))
	mm.LocalEnergyMinimizer.minimize(context)
	sampler_state.update_from_context(context)
	final_energy = thermodynamic_state.reduced_potential(context)
	print("Sampler state {}: final energy {:8.3f}kT".format(k, final_energy))
	del context
print('Minimized......')

if runflag == 'run':
	repex_simulation = ReplicaExchangeSampler(mcmc_moves=lsd_move, number_of_iterations=productionIterations)

tmp_dir = './trajectory/'
storage = os.path.join(tmp_dir, 'solvent.nc')
reporter = MultiStateReporter(storage, checkpoint_interval=iterationsPerCheckpoint)
if runflag != 'run':
    repex_simulation = ReplicaExchangeSampler.from_storage(reporter)
else:
    repex_simulation.create(thermodynamic_states, sampler_states, reporter)
    print('Equilibrating......')
    repex_simulation.equilibrate(equilibrationIterations)
    print('Simulating......')
if runflag == 'recover' or runflag == 'run':
        repex_simulation.run()
elif runflag == 'extend':
        repex_simulation.extend(extendIterations)

#will add all iterations even if coming from a previous restart
all_iters = repex_simulation.iteration
print('All iterations = {}'.format(all_iters))

analyzer = ReplicaExchangeAnalyzer(reporter)
iterations_to_analyze=int(all_iters/10)
print('Iterations to analyze = {}'.format(iterations_to_analyze))
for i in range(1, iterations_to_analyze+1):
	samples_to_analyze=i*10
	analyzer.max_n_iterations = samples_to_analyze
	Delta_f_ij, dDelta_f_ij = analyzer.get_free_energy()
	print("Relative free energy change in {0} = {1} +- {2}"
		.format('solvent', Delta_f_ij[0, nstates - 1]*kTtokcal, dDelta_f_ij[0, nstates - 1]*kTtokcal))

[matrix,eigenvalues,ineff]=analyzer.generate_mixing_statistics()
print("Mixing Stats")
print(matrix)
print(eigenvalues)
print(ineff)
