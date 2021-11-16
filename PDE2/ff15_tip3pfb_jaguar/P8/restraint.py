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
waterModel='tip3pfb'

#System
ionicStrength=0.15*unit.molar
boxPadding=0.5*unit.nanometers

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

pdb = app.PDBFile('input/complex.pdb')
forcefield = app.ForceField(proteinXML, waterXML, ionXML, forcefieldXML, 'input/ligand.xml')
solvated = app.Modeller (pdb.topology, pdb.positions)
solvated.addSolvent(forcefield, model='tip3p', ionicStrength=ionicStrength, neutralize=True, padding=boxPadding)
system = forcefield.createSystem(solvated.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005, hydrogenMass=4*unit.amus)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))
forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
nforces = system.getNumForces()

restraint_state = RestraintComposableState(lambda_restraints=1.0)
ligand_a_atoms = [5332,5333,5334,5335,5336,5337,5338,5339,5340,5341,5342,5343,5344,5345,5346,5347,5348,5349,5350,5351,5352,5353,5354,5355,5356,5357,5358,5359,5360,5361,5362,5363,5364,5365,5366,5367,5368,5369,5370,5371,5372,5373,5374,5375,5376,5377,5378,5379,5380]
ligand_a_bonds = []
ligand_a_angles =[]
ligand_a_torsions =[]

#Apply harmonic restraints
springconstant=1*unit.kilocalorie_per_mole/unit.angstroms**2
bondlength=2.42591*unit.angstroms
bond_energy_function = "lambda_restraints*(K/2)*(r-r0)^2;"
harmonicforce=mm.CustomBondForce(bond_energy_function)
harmonicforce.addPerBondParameter('r0')
harmonicforce.addPerBondParameter('K')
harmonicforce.addGlobalParameter('lambda_restraints', 1.0)
harmonicforce.addBond(3548,5357, [bondlength,springconstant])
new_force_index=system.getNumForces()
harmonicforce.setForceGroup(new_force_index)
system.addForce(harmonicforce)

angleforceconstant=5*unit.kilocalorie_per_mole/unit.radian**2
restrainedangle=1.58074*unit.radian
angle_energy_function = "lambda_restraints*0.5*k*(theta-theta0)^2"
angleforce=mm.CustomAngleForce(angle_energy_function)
angleforce.addPerAngleParameter('theta0')
angleforce.addPerAngleParameter('k')
angleforce.addGlobalParameter('lambda_restraints', 1.0)
angleforce.addAngle(5341,5357,3548, [restrainedangle,angleforceconstant])
angleforce.setForceGroup(new_force_index)

restrainedangle=1.93988*unit.radian
angleforce.addAngle(5357,3548,3545, [restrainedangle,angleforceconstant])
system.addForce(angleforce)

torsionforceconstant=5*unit.kilocalorie_per_mole
torsion_energy_function = "lambda_restraints*k*(1+cos(n*theta-theta0))"
torsionforce=mm.CustomTorsionForce(torsion_energy_function)
torsionforce.addPerTorsionParameter('n')
torsionforce.addPerTorsionParameter('theta0')
torsionforce.addPerTorsionParameter('k')
torsionforce.addGlobalParameter('lambda_restraints', 1.0)
torsionforce.addTorsion(5349,5341,5357,3548, [1,-0.0944614*unit.radian,torsionforceconstant])
torsionforce.setForceGroup(new_force_index)

torsionforce.addTorsion(5341,5357,3548,3545, [1,-4.30348*unit.radian,torsionforceconstant])

torsionforce.addTorsion(5357,3548,3545,3547, [1,-5.69933*unit.radian,torsionforceconstant])
system.addForce(torsionforce)  
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
#file = open('DEBUG_restraint.xml','w')
#file.write(mm.XmlSerializer.serialize(sys))
#file.close()

nstates=restraintSteps+1
print("There will be ", nstates, " states in total")
print("Lambda exponent: ", restraintExponent)
print("stepsPerIteration:", stepsPerIteration, " productionIterations: ", productionIterations, "equilibrationIterations: ", equilibrationIterations)
print("Timestep: ", timeStep)
box_vec = alchemical_system_in.getDefaultPeriodicBoxVectors()
print("Box vectors:", box_vec)

#Setup lambda schedule
lambdas = np.linspace(1.0, 0.0, nstates)

for j in range(restraintSteps):
        if(lambdas[j] < restraintInversion):
                lambdas[j]=restraintInversion*pow(lambdas[j]/restraintInversion,1/restraintExponent)
        elif(lambdas[j] > restraintInversion):
                lambdas[j]=restraintInversion+(1-restraintInversion)*pow((lambdas[j]-restraintInversion)/(1-restraintInversion),restraintExponent)

#Sanity check
print("")
print("Lambdas matrix")
print("Lrestraints")
for j in range(len(lambdas)):
    print("%-15.2f" % lambdas[j])

sampler_states = list()
thermodynamic_states = list()

for k in range(nstates):
    compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
    compound_state.lambda_sterics_A=1.0
    compound_state.lambda_electrostatics_A=1.0
    compound_state.lambda_restraints=lambdas[k]
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
storage = os.path.join(tmp_dir, 'restraint.nc')
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
              .format('restraint', Delta_f_ij[0, nstates - 1]*kTtokcal, dDelta_f_ij[0, nstates - 1]*kTtokcal))

[matrix,eigenvalues,ineff]=analyzer.generate_mixing_statistics()
print("Mixing Stats")
print(matrix)
print(eigenvalues)
print(ineff)
