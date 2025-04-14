# OpenMM Python API Specification

This document describes the overview of OpenMM's Python API, its main classes, usage, and use cases. It also includes details of arguments for each function and method.

## Table of Contents

1. [Introduction](#introduction)
2. [Module Structure](#module-structure)
3. [Main Components](#main-components)
   - [System](#system)
   - [Force](#force)
   - [Integrator](#integrator)
   - [Context](#context)
   - [Platform](#platform)
   - [Modeller](#modeller-openmmappmodeller)
4. [Setting Up Simulations](#setting-up-simulations)
5. [File Input and Output](#file-input-and-output)
6. [Analysis and Visualization](#analysis-and-visualization)
7. [Advanced Features](#advanced-features)
8. [Code Examples](#code-examples)
9. [Troubleshooting](#troubleshooting)
10. [Function Reference](#function-reference)

## Introduction

OpenMM is a toolkit for molecular simulation. It can be used as a standalone application, or as a library called from your own code. It provides a combination of extreme flexibility, openness, and high performance (especially on recent GPUs) through its custom forces and integrators, making it truly unique among simulation codes.

## Module Structure

The Python API of OpenMM consists of the following main modules:

- `openmm` - Interface to the core library
- `openmm.app` - Tools for building molecular systems, file I/O, and running simulations
- `openmm.unit` - Working with unit-bearing physical quantities

## Main Components

### System

The `System` class represents a physical description of a molecular system. It consists of particles (atoms), forces, and constraints.

#### Main Methods

**addParticle**
```python
System.addParticle(mass)
```
- `mass` (quantified mass) - The mass of the particle (with units)

**addConstraint**
```python
System.addConstraint(particle1, particle2, distance)
```
- `particle1` (int) - Index of the first particle
- `particle2` (int) - Index of the second particle
- `distance` (quantified length) - The constrained distance (with units)

**addForce**
```python
System.addForce(force)
```
- `force` (Force) - Force object to add to the system

**getNumParticles**
```python
System.getNumParticles()
```
- Returns: Number of particles in the system

**getNumConstraints**
```python
System.getNumConstraints()
```
- Returns: Number of constraints in the system

**getNumForces**
```python
System.getNumForces()
```
- Returns: Number of forces in the system

**getForce**
```python
System.getForce(index)
```
- `index` (int) - Index of the force
- Returns: Force object with the specified index

**setDefaultPeriodicBoxVectors**
```python
System.setDefaultPeriodicBoxVectors(a, b, c)
```
- `a`, `b`, `c` (Vec3) - Periodic boundary condition box vectors (with units)

### Force

`Force` is an abstract base class representing interactions in a molecular system. OpenMM has numerous built-in force classes:

#### HarmonicBondForce

```python
force = HarmonicBondForce()
```

**addBond**
```python
HarmonicBondForce.addBond(particle1, particle2, length, k)
```
- `particle1` (int) - Index of the first particle
- `particle2` (int) - Index of the second particle
- `length` (quantified length) - Equilibrium bond length (with units)
- `k` (quantified energy/length^2) - Force constant (with units)

#### HarmonicAngleForce

```python
force = HarmonicAngleForce()
```

**addAngle**
```python
HarmonicAngleForce.addAngle(particle1, particle2, particle3, angle, k)
```
- `particle1` (int) - Index of the first particle
- `particle2` (int) - Index of the central particle
- `particle3` (int) - Index of the third particle
- `angle` (quantified angle) - Equilibrium angle (with units, radians)
- `k` (quantified energy/angle^2) - Force constant (with units)

#### PeriodicTorsionForce

```python
force = PeriodicTorsionForce()
```

**addTorsion**
```python
PeriodicTorsionForce.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
```
- `particle1`, `particle2`, `particle3`, `particle4` (int) - Indices of the four particles
- `periodicity` (int) - Periodicity
- `phase` (quantified angle) - Phase angle (with units, radians)
- `k` (quantified energy) - Force constant (with units)

#### NonbondedForce

```python
force = NonbondedForce()
```

**addParticle**
```python
NonbondedForce.addParticle(charge, sigma, epsilon)
```
- `charge` (quantified charge) - Particle charge (with units)
- `sigma` (quantified length) - Lennard-Jones sigma parameter (with units)
- `epsilon` (quantified energy) - Lennard-Jones epsilon parameter (with units)

**addException**
```python
NonbondedForce.addException(particle1, particle2, chargeProd, sigma, epsilon, replace=False)
```
- `particle1`, `particle2` (int) - Indices of two particles
- `chargeProd` (quantified charge^2) - Product of charges (with units)
- `sigma` (quantified length) - Lennard-Jones sigma parameter (with units)
- `epsilon` (quantified energy) - Lennard-Jones epsilon parameter (with units)
- `replace` (bool) - Whether to replace existing exceptions

**setNonbondedMethod**
```python
NonbondedForce.setNonbondedMethod(method)
```
- `method` (int) - Method for calculating nonbonded interactions (`NoCutoff`, `CutoffNonPeriodic`, `CutoffPeriodic`, `Ewald`, `PME`, `LJPME`)

**setCutoffDistance**
```python
NonbondedForce.setCutoffDistance(distance)
```
- `distance` (quantified length) - Cutoff distance (with units)

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

**addPerBondParameter**
```python
CustomBondForce.addPerBondParameter(name)
```
- `name` (str) - Parameter name

**addBond**
```python
CustomBondForce.addBond(particle1, particle2, parameters)
```
- `particle1`, `particle2` (int) - Indices of two particles
- `parameters` (list) - List of parameter values

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

**addPerParticleParameter**
```python
CustomNonbondedForce.addPerParticleParameter(name)
```
- `name` (str) - Parameter name

**addParticle**
```python
CustomNonbondedForce.addParticle(parameters)
```
- `parameters` (list) - List of parameter values

**addExclusion**
```python
CustomNonbondedForce.addExclusion(particle1, particle2)
```
- `particle1`, `particle2` (int) - Indices of two particles to exclude

### Integrator

The `Integrator` class implements numerical integration of the equations of motion. Main integrators include:

#### VerletIntegrator

```python
integrator = VerletIntegrator(stepSize)
```
- `stepSize` (quantified time) - Integration time step (with units)

#### LangevinIntegrator

```python
integrator = LangevinIntegrator(temperature, frictionCoeff, stepSize)
```
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frictionCoeff` (quantified inverse time) - Friction coefficient (with units)
- `stepSize` (quantified time) - Integration time step (with units)

#### BrownianIntegrator

```python
integrator = BrownianIntegrator(temperature, frictionCoeff, stepSize)
```
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frictionCoeff` (quantified inverse time) - Friction coefficient (with units)
- `stepSize` (quantified time) - Integration time step (with units)

#### VariableVerletIntegrator

```python
integrator = VariableVerletIntegrator(errorTol)
```
- `errorTol` (float) - Error tolerance (dimensionless)

#### VariableLangevinIntegrator

```python
integrator = VariableLangevinIntegrator(temperature, frictionCoeff, errorTol)
```
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frictionCoeff` (quantified inverse time) - Friction coefficient (with units)
- `errorTol` (float) - Error tolerance (dimensionless)

#### MTSIntegrator

```python
integrator = MTSIntegrator(timestep, loops)
```
- `timestep` (list) - List of hierarchical time steps (with units)
- `loops` (list) - List of the number of loops to perform at each level

#### AMDIntegrator

```python
integrator = AMDIntegrator(temperature, frictionCoeff, stepSize, alpha, E)
```
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frictionCoeff` (quantified inverse time) - Friction coefficient (with units)
- `stepSize` (quantified time) - Integration time step (with units)
- `alpha` (quantified energy) - Boost parameter (with units)
- `E` (quantified energy) - Energy threshold (with units)

### Context

`Context` holds the runtime state of a simulation. It is created by combining a specific System, Integrator, and Platform.

```python
context = Context(system, integrator, platform=None, properties=None)
```
- `system` (System) - System to simulate
- `integrator` (Integrator) - Integrator to use
- `platform` (Platform, optional) - Platform to use
- `properties` (dict, optional) - Platform-specific properties

#### Main Methods

**getState**
```python
Context.getState(getPositions=False, getVelocities=False, getForces=False, 
                getEnergy=False, getParameters=False, getParameterDerivatives=False, 
                getIntegratorParameters=False, enforcePeriodicBox=False, groups=-1)
```
- `getPositions` (bool) - Whether to include positions
- `getVelocities` (bool) - Whether to include velocities
- `getForces` (bool) - Whether to include forces
- `getEnergy` (bool) - Whether to include energy
- `getParameters` (bool) - Whether to include context parameters
- `getParameterDerivatives` (bool) - Whether to include parameter derivatives
- `getIntegratorParameters` (bool) - Whether to include integrator parameters
- `enforcePeriodicBox` (bool) - Whether to force positions within the periodic box
- `groups` (int) - Force groups to include (bitmask)

**setState**
```python
Context.setState(state)
```
- `state` (State) - State to set

**reinitialize**
```python
Context.reinitialize(preserveState=False)
```
- `preserveState` (bool) - Whether to preserve the current state

**setPositions**
```python
Context.setPositions(positions)
```
- `positions` (list) - List of particle positions (with units)

**setVelocities**
```python
Context.setVelocities(velocities)
```
- `velocities` (list) - List of particle velocities (with units)

**setParameter**
```python
Context.setParameter(name, value)
```
- `name` (str) - Parameter name
- `value` (float) - Parameter value

**getParameter**
```python
Context.getParameter(name)
```
- `name` (str) - Parameter name
- Returns: Parameter value

**setPeriodicBoxVectors**
```python
Context.setPeriodicBoxVectors(a, b, c)
```
- `a`, `b`, `c` (Vec3) - Periodic boundary condition box vectors (with units)

### Platform

The `Platform` class represents a hardware platform (CPU, CUDA, OpenCL, Reference) to run simulations on.

#### Main Methods

**getName**
```python
Platform.getName()
```
- Returns: Platform name

**getPropertyNames**
```python
Platform.getPropertyNames()
```
- Returns: List of property names defined for this platform

**getPropertyValue**
```python
Platform.getPropertyValue(context, property)
```
- `context` (Context) - Context
- `property` (str) - Property name
- Returns: Value of the property

**setPropertyValue**
```python
Platform.setPropertyValue(context, property, value)
```
- `context` (Context) - Context
- `property` (str) - Property name
- `value` (str) - Value to set

#### Static Methods

**getPlatformByName**
```python
Platform.getPlatformByName(name)
```
- `name` (str) - Platform name
- Returns: Named platform

**getNumPlatforms**
```python
Platform.getNumPlatforms()
```
- Returns: Number of available platforms

**findPlatform**
```python
Platform.findPlatform(capabilities)
```
- `capabilities` (list) - List of required capabilities
- Returns: Platform meeting the criteria

### Modeller (openmm.app.Modeller)

The `Modeller` class provides tools for editing molecular models, enabling operations such as adding water and hydrogen atoms.

```python
modeller = Modeller(topology, positions)
```
- `topology` (Topology) - Topology of the initial molecular system
- `positions` (list) - Initial atomic positions (with units)

#### Basic Properties and Methods

**topology**
```python
modeller.topology
```
- Topology object describing the structure of the system

**positions**
```python
modeller.positions
```
- List of atomic positions (with units)

**getTopology**
```python
modeller.getTopology()
```
- Returns: The model's Topology

**getPositions**
```python
modeller.getPositions()
```
- Returns: List of atomic positions (with units)

#### Editing Molecular Structure

**add**
```python
modeller.add(addTopology, addPositions)
```
- `addTopology` (Topology) - Topology of the molecular system to add
- `addPositions` (list) - Atomic positions of the molecular system to add (with units)
- Description: Adds new molecular structures to the existing model

**delete**
```python
modeller.delete(toDelete)
```
- `toDelete` (list) - List of atoms, residues, or chains to delete
- Description: Removes specified atoms, residues, or chains from the model

**deleteWater**
```python
modeller.deleteWater()
```
- Description: Removes all water molecules from the model

**convertWater**
```python
modeller.convertWater(model='tip3p')
```
- `model` (str) - Target water model ('tip3p', 'tip4pew', 'tip5p', 'spce', 'swm4ndp', etc.)
- Description: Converts water molecules to the specified model. Supports conversion between models with different virtual site structures

#### Solvation and Ion Addition

**addSolvent**
```python
modeller.addSolvent(forcefield, model='tip3p', boxSize=None, boxVectors=None, 
                   padding=None, numAdded=None, boxShape='cube', 
                   positiveIon='Na+', negativeIon='Cl-', 
                   ionicStrength=0*molar, neutralize=True, residueTemplates=dict())
```
- `forcefield` (ForceField) - Force field to use
- `model` (str) - Water model to use ('tip3p', 'spce', 'tip4pew', 'tip5p', 'swm4ndp')
- `boxSize` (Vec3, optional) - Size of rectangular box (with units)
- `boxVectors` (tuple, optional) - Periodic boundary condition box vectors
- `padding` (quantity, optional) - Padding distance of water around the solute (with units)
- `numAdded` (int, optional) - Total number of solvent molecules (water + ions) to add
- `boxShape` (str) - Box shape ('cube', 'dodecahedron', 'octahedron')
- `positiveIon` (str) - Type of positive ion ('Cs+', 'K+', 'Li+', 'Na+', 'Rb+')
- `negativeIon` (str) - Type of negative ion ('Cl-', 'Br-', 'F-', 'I-')
- `ionicStrength` (quantity) - Ionic strength (with units, molar)
- `neutralize` (bool) - Whether to neutralize the system
- `residueTemplates` (dict) - Specification of templates to use for specific residues
- Description: Solvates the system by adding water molecules and ions. Supports various box shapes and water models

#### Adding Hydrogen Atoms

**addHydrogens**
```python
modeller.addHydrogens(forcefield=None, pH=7.0, variants=None, platform=None, residueTemplates=dict())
```
- `forcefield` (ForceField, optional) - Force field to use
- `pH` (float) - pH value (used to determine protonation states)
- `variants` (list, optional) - List of residue variants (HIE, HID, HIP, etc.)
- `platform` (Platform, optional) - Platform to use for calculating hydrogen atom positions
- `residueTemplates` (dict) - Specification of templates to use for specific residues
- Description: Adds hydrogen atoms to the model. Can automatically select protonation states based on pH
- Returns: List of variants selected for each residue

#### Special Particle Operations

**addExtraParticles**
```python
modeller.addExtraParticles(forcefield, ignoreExternalBonds=False, residueTemplates=dict())
```
- `forcefield` (ForceField) - Force field to use
- `ignoreExternalBonds` (bool) - Whether to ignore external bonds
- `residueTemplates` (dict) - Specification of templates to use for specific residues
- Description: Adds additional particles to the model, such as virtual sites defined in the force field

#### Building Membrane Systems

**addMembrane**
```python
modeller.addMembrane(forcefield, lipidType='POPC', membraneCenterZ=0*nanometer, 
                    minimumPadding=1*nanometer, positiveIon='Na+', negativeIon='Cl-',
                    ionicStrength=0*molar, neutralize=True, residueTemplates=dict(), platform=None)
```
- `forcefield` (ForceField) - Force field to use
- `lipidType` (str) - Type of lipid to use (e.g., 'POPC')
- `membraneCenterZ` (quantity) - Z-coordinate of the membrane center (with units)
- `minimumPadding` (quantity) - Minimum distance between membrane edge and solute molecules (with units)
- `positiveIon`, `negativeIon` (str) - Types of ions to use
- `ionicStrength` (quantity) - Ionic strength (with units, molar)
- `neutralize` (bool) - Whether to neutralize the system
- `residueTemplates` (dict) - Specification of templates to use for specific residues
- `platform` (Platform, optional) - Platform to use for energy minimization
- Description: Adds a lipid bilayer to the molecular system. Useful for building protein-membrane complexes

#### Helper Methods

**loadHydrogenDefinitions**
```python
Modeller.loadHydrogenDefinitions(file)
```
- `file` (str) - Path to hydrogen definition XML file
- Description: Loads custom hydrogen definitions (class method)

### Force

`Force` is an abstract base class representing interactions in a molecular system. OpenMM has numerous built-in force classes:

#### HarmonicBondForce

```python
force = HarmonicBondForce()
```

**addBond**
```python
HarmonicBondForce.addBond(particle1, particle2, length, k)
```
- `particle1` (int) - Index of the first particle
- `particle2` (int) - Index of the second particle
- `length` (quantified length) - Equilibrium bond length (with units)
- `k` (quantified energy/length^2) - Force constant (with units)

#### HarmonicAngleForce

```python
force = HarmonicAngleForce()
```

**addAngle**
```python
HarmonicAngleForce.addAngle(particle1, particle2, particle3, angle, k)
```
- `particle1` (int) - Index of the first particle
- `particle2` (int) - Index of the central particle
- `particle3` (int) - Index of the third particle
- `angle` (quantified angle) - Equilibrium angle (with units, radians)
- `k` (quantified energy/angle^2) - Force constant (with units)

#### PeriodicTorsionForce

```python
force = PeriodicTorsionForce()
```

**addTorsion**
```python
PeriodicTorsionForce.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
```
- `particle1`, `particle2`, `particle3`, `particle4` (int) - Indices of the four particles
- `periodicity` (int) - Periodicity
- `phase` (quantified angle) - Phase angle (with units, radians)
- `k` (quantified energy) - Force constant (with units)

#### NonbondedForce

```python
force = NonbondedForce()
```

**addParticle**
```python
NonbondedForce.addParticle(charge, sigma, epsilon)
```
- `charge` (quantified charge) - Particle charge (with units)
- `sigma` (quantified length) - Lennard-Jones sigma parameter (with units)
- `epsilon` (quantified energy) - Lennard-Jones epsilon parameter (with units)

**addException**
```python
NonbondedForce.addException(particle1, particle2, chargeProd, sigma, epsilon, replace=False)
```
- `particle1`, `particle2` (int) - Indices of two particles
- `chargeProd` (quantified charge^2) - Product of charges (with units)
- `sigma` (quantified length) - Lennard-Jones sigma parameter (with units)
- `epsilon` (quantified energy) - Lennard-Jones epsilon parameter (with units)
- `replace` (bool) - Whether to replace existing exceptions

**setNonbondedMethod**
```python
NonbondedForce.setNonbondedMethod(method)
```
- `method` (int) - Method for calculating nonbonded interactions (`NoCutoff`, `CutoffNonPeriodic`, `CutoffPeriodic`, `Ewald`, `PME`, `LJPME`)

**setCutoffDistance**
```python
NonbondedForce.setCutoffDistance(distance)
```
- `distance` (quantified length) - Cutoff distance (with units)

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

**addPerBondParameter**
```python
CustomBondForce.addPerBondParameter(name)
```
- `name` (str) - Parameter name

**addBond**
```python
CustomBondForce.addBond(particle1, particle2, parameters)
```
- `particle1`, `particle2` (int) - Indices of two particles
- `parameters` (list) - List of parameter values

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

**addPerParticleParameter**
```python
CustomNonbondedForce.addPerParticleParameter(name)
```
- `name` (str) - Parameter name

**addParticle**
```python
CustomNonbondedForce.addParticle(parameters)
```
- `parameters` (list) - List of parameter values

**addExclusion**
```python
CustomNonbondedForce.addExclusion(particle1, particle2)
```
- `particle1`, `particle2` (int) - Indices of two particles to exclude

## Setting Up Simulations

The most common way to set up a simulation in OpenMM is to use the `openmm.app.Simulation` class.

```python
simulation = Simulation(topology, system, integrator, platform=None, platformProperties=None, state=None)
```
- `topology` (Topology) - Topology of the molecular system
- `system` (System) - System to simulate
- `integrator` (Integrator) - Integrator to use
- `platform` (Platform, optional) - Platform to use
- `platformProperties` (dict, optional) - Platform-specific properties
- `state` (State, optional) - Initial state

### Main Methods

**context**
```python
Simulation.context
```
- Context object for the simulation

**minimizeEnergy**
```python
Simulation.minimizeEnergy(tolerance=10*kilojoules/mole, maxIterations=0)
```
- `tolerance` (quantified energy) - Convergence tolerance (with units)
- `maxIterations` (int) - Maximum number of iterations (0 means unlimited)

**step**
```python
Simulation.step(steps)
```
- `steps` (int) - Number of steps to run

**saveState**
```python
Simulation.saveState(file)
```
- `file` (str or file-like object) - File to save the state to

**loadState**
```python
Simulation.loadState(file)
```
- `file` (str or file-like object) - File to load the state from

**reporters**
```python
Simulation.reporters
```
- List of reporters

The basic setup for a simulation is done in the following steps:

1. Define topology (`Topology`) and position information
2. Select a force field (`ForceField`)
3. Create a system (`System`)
4. Choose an integrator (`Integrator`)
5. Create a simulation (`Simulation`) object
6. Set positions and velocities
7. Add reporters
8. Run the simulation

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# 1. Load topology and position information from a PDB file
pdb = PDBFile('input.pdb')

# 2. Specify force field files
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# 3. Create a system
system = forcefield.createSystem(
    pdb.topology, 
    nonbondedMethod=PME,  # Use periodic boundary conditions
    nonbondedCutoff=1*nanometer,  # Cutoff for nonbonded interactions
    constraints=HBonds    # Constrain hydrogen bond lengths
)

# 4. Choose an integrator (here Langevin dynamics)
integrator = LangevinIntegrator(
    300*kelvin,       # Temperature
    1/picosecond,     # Friction coefficient
    0.002*picoseconds # Time step
)

# 5. Create the simulation
simulation = Simulation(pdb.topology, system, integrator)

# 6. Set initial positions
simulation.context.setPositions(pdb.positions)

# 7. Add reporters
simulation.reporters.append(PDBReporter('output.pdb', 1000))  # Output structure every 1000 steps
simulation.reporters.append(StateDataReporter(
    'output.csv', 1000, step=True, temperature=True, potentialEnergy=True))

# 8. Perform energy minimization
simulation.minimizeEnergy()

# 9. Run the simulation (e.g., 100ps of molecular dynamics)
simulation.step(50000)  # 0.002 ps × 50,000 steps = 100 ps
```

## File Input and Output

OpenMM supports reading and writing various molecular file formats:

### Input Files

#### PDBFile

```python
pdb = PDBFile(file)
```
- `file` (str or file-like object) - Path or object of the PDB file

#### AmberPrmtopFile / AmberInpcrdFile

```python
prmtop = AmberPrmtopFile(file)
inpcrd = AmberInpcrdFile(file)
```
- `file` (str or file-like object) - Path or object of the AMBER file

#### GromacsTopFile / GromacsGroFile

```python
top = GromacsTopFile(file, periodicBoxVectors=None, includeDir=None)
gro = GromacsGroFile(file)
```
- `file` (str or file-like object) - Path or object of the GROMACS file
- `periodicBoxVectors` (tuple, optional) - Periodic boundary condition box vectors
- `includeDir` (str, optional) - Directory for include files

#### CharmmPsfFile / CharmmCrdFile

```python
psf = CharmmPsfFile(file)
crd = CharmmCrdFile(file)
```
- `file` (str or file-like object) - Path or object of the CHARMM file

### Output/Reporters

#### PDBReporter

```python
reporter = PDBReporter(file, reportInterval, enforcePeriodicBox=True)
```
- `file` (str or file-like object) - Path or object of the output file
- `reportInterval` (int) - Reporting interval (number of steps)
- `enforcePeriodicBox` (bool) - Whether to fit particles within the periodic box

#### DCDReporter

```python
reporter = DCDReporter(file, reportInterval, append=False, enforcePeriodicBox=True)
```
- `file` (str or file-like object) - Path or object of the output file
- `reportInterval` (int) - Reporting interval (number of steps)
- `append` (bool) - Whether to append to existing file
- `enforcePeriodicBox` (bool) - Whether to fit particles within the periodic box

#### XTCReporter

```python
reporter = XTCReporter(file, reportInterval, append=False, enforcePeriodicBox=True)
```
- `file` (str or file-like object) - Path or object of the output file
- `reportInterval` (int) - Reporting interval (number of steps)
- `append` (bool) - Whether to append to existing file
- `enforcePeriodicBox` (bool) - Whether to fit particles within the periodic box

#### StateDataReporter

```python
reporter = StateDataReporter(file, reportInterval, step=False, time=False, 
                            potentialEnergy=False, kineticEnergy=False, 
                            totalEnergy=False, temperature=False, volume=False, 
                            density=False, progress=False, remainingTime=False, 
                            speed=False, elapsedTime=False, separator=',', 
                            systemMass=None, totalSteps=None)
```
- `file` (str, file-like object, or None) - Output file or console output (None)
- `reportInterval` (int) - Reporting interval (number of steps)
- Other arguments - Specify physical quantities and output formats to report

#### CheckpointReporter

```python
reporter = CheckpointReporter(file, reportInterval)
```
- `file` (str) - Path of the output file
- `reportInterval` (int) - Reporting interval (number of steps)

## Analysis and Visualization

OpenMM provides several tools for analyzing simulation results:

- Use `StateDataReporter` to record time series of physical quantities such as energy and temperature
- Output trajectories in PDB, DCD, XTC formats for analysis with external tools like VMD, PyMOL, MDAnalysis
- Directly retrieve information such as positions, velocities, and forces from the simulation's `State` object for analysis

### State Object

The `State` object represents an instantaneous state of a simulation.

```python
state = context.getState(getPositions=False, getVelocities=False, getForces=False, 
                        getEnergy=False, getParameters=False, getParameterDerivatives=False, 
                        enforcePeriodicBox=False, groups=-1)
```

#### Main Methods

**getPositions**
```python
State.getPositions(asNumpy=False)
```
- `asNumpy` (bool) - Whether to return as a NumPy array
- Returns: Particle positions (with units)

**getVelocities**
```python
State.getVelocities(asNumpy=False)
```
- `asNumpy` (bool) - Whether to return as a NumPy array
- Returns: Particle velocities (with units)

**getForces**
```python
State.getForces(asNumpy=False)
```
- `asNumpy` (bool) - Whether to return as a NumPy array
- Returns: Forces acting on particles (with units)

**getPotentialEnergy**
```python
State.getPotentialEnergy()
```
- Returns: Potential energy (with units)

**getKineticEnergy**
```python
State.getKineticEnergy()
```
- Returns: Kinetic energy (with units)

**getPeriodicBoxVectors**
```python
State.getPeriodicBoxVectors()
```
- Returns: Periodic boundary condition box vectors (with units)

## Advanced Features

### Custom Forces

One of the most powerful features of OpenMM is the ability to define custom forces. These allow you to express almost any physical interaction:

#### CustomBondForce

```python
force = CustomBondForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom bond forces define a potential with an arbitrary functional form between two atoms. The expression can use the variable `r` (bond length) and user-defined parameters.

#### CustomAngleForce

```python
force = CustomAngleForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom angle forces define a potential with an arbitrary functional form among three atoms. The expression can use the variable `theta` (angle) and user-defined parameters.

#### CustomTorsionForce

```python
force = CustomTorsionForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom torsion forces define a potential with an arbitrary functional form among four atoms. The expression can use the variable `theta` (dihedral angle) and user-defined parameters.

#### CustomNonbondedForce

```python
force = CustomNonbondedForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom nonbonded forces define a potential with an arbitrary functional form between particle pairs. The expression can use the variable `r` (distance between particles) and user-defined parameters.

#### CustomExternalForce

```python
force = CustomExternalForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom external forces define an external potential with an arbitrary functional form acting on individual particles. The expression can use the variables `x`, `y`, `z` (particle coordinates) and user-defined parameters.

### Constrained Simulations

OpenMM includes tools for applying various types of constraints:

#### Bond Length Constraints

Bond length constraints are specified through the `constraints` parameter of the `System.createSystem()` method:

```python
system = forcefield.createSystem(topology, constraints=HBonds)
```

- `constraints` - Constraint level (`None`, `HBonds`, `AllBonds`, `HAngles`)

#### CustomCVForce

```python
force = CustomCVForce(energy)
```
- `energy` (str) - Mathematical expression defining the potential energy

Custom collective variable forces define a potential based on complex collective variables (CVs).

#### RMSDForce

```python
force = RMSDForce(referencePositions, particles=None)
```
- `referencePositions` (list) - List of reference positions (with units)
- `particles` (list, optional) - List of particle indices to include

RMSD-based constraint forces define a potential based on the RMSD (root mean square deviation) between the current structure and a reference structure.

### Advanced Sampling Methods

OpenMM includes various techniques for improving sampling of systems that cannot be efficiently explored with regular molecular dynamics simulations. Below are the main advanced sampling methods and examples of their use.

#### MTSIntegrator (Multiple Time Step Integration)

```python
integrator = MTSIntegrator(timestep, loops)
```
- `timestep` (list) - List of hierarchical time steps (with units)
- `loops` (list) - List of the number of loops to perform at each level

The multiple time step integrator improves computational efficiency by evaluating different forces at different frequencies. For example, it distinguishes between rapidly changing forces (like bond forces) and slowly changing forces (like long-range interactions).

**Usage Example**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# Load structure from PDB file
pdb = PDBFile('input.pdb')

# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system and divide forces into two groups
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Find NonbondedForce and set it to group 1
for i in range(system.getNumForces()):
    force = system.getForce(i)
    if isinstance(force, NonbondedForce):
        force.setForceGroup(1)
    else:
        # Set other forces to group 0
        force.setForceGroup(0)

# Set up MTSIntegrator
# Group 0 (bond forces, etc.) evaluated at 0.5 fs
# Group 1 (nonbonded forces, etc.) evaluated at 2.0 fs
inner_ts = 0.5*femtoseconds
outer_ts = 2.0*femtoseconds

integrator = MTSIntegrator([inner_ts, outer_ts], [4, 1])
# Group 0 is evaluated 4 times with inner_ts, and Group 1 is evaluated once with outer_ts

# Prepare simulation
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Energy minimization
simulation.minimizeEnergy()

# Run simulation
simulation.reporters.append(StateDataReporter('mts_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))
simulation.step(50000)  # Run 50,000 steps of simulation
```

#### AMDIntegrator (Accelerated Molecular Dynamics)

```python
integrator = AMDIntegrator(temperature, frictionCoeff, stepSize, alpha, E)
```
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frictionCoeff` (quantified inverse time) - Friction coefficient (with units)
- `stepSize` (quantified time) - Integration time step (with units)
- `alpha` (quantified energy) - Boost parameter (with units)
- `E` (quantified energy) - Energy threshold (with units)

The accelerated molecular dynamics (AMD) integrator accelerates exploration of a molecular system's conformational space by modifying the potential energy landscape. It is particularly effective for systems with high energy barriers where transitions are unlikely in regular simulations.

**Usage Example**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np

# Load structure from PDB file
pdb = PDBFile('protein.pdb')

# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# First run a short standard simulation to estimate energy
temp_integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
temp_simulation = Simulation(pdb.topology, system, temp_integrator)
temp_simulation.context.setPositions(pdb.positions)
temp_simulation.minimizeEnergy()
temp_simulation.context.setVelocitiesToTemperature(300*kelvin)
temp_simulation.step(5000)  # Run a short simulation

# Calculate average potential energy
state = temp_simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

# Set AMD parameters
# E_thresh is set slightly higher than the system's average potential energy
E_thresh = potential_energy * 1.05 * kilojoule_per_mole
# alpha is set to about 1/5 of E_thresh
alpha = 0.2 * E_thresh

# Use AMDIntegrator
amd_integrator = AMDIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, alpha, E_thresh)

# Prepare new simulation
simulation = Simulation(pdb.topology, system, amd_integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

# Add reporters
simulation.reporters.append(DCDReporter('amd_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('amd_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# Run simulation
simulation.step(500000)  # Run 500,000 steps of simulation (1 ns)
```

#### MetaDynamics

```python
meta = MetaDynamics(system, variables, temperature, biasFactor, height, frequency)
```
- `system` (System) - System to simulate
- `variables` (list) - List of collective variables
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `biasFactor` (float) - Bias factor
- `height` (quantified energy) - Energy height of Gaussian hills (with units)
- `frequency` (int) - Frequency of adding hills (number of steps)

Metadynamics improves sampling of molecular systems by applying energy penalties to regions of conformational space that have already been visited. This is particularly useful for studying complex structural changes such as protein folding or conformational transitions.

**Usage Example**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# Load structure from PDB file
pdb = PDBFile('protein.pdb')

# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Create basic integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Define a specific atom distance (collective variable)
# Example: Distance between CA atoms of residues 10 and 20
atoms = list(pdb.topology.atoms())
ca_atoms = [atom.index for atom in atoms if atom.name == 'CA']

# Find CA atoms of residues 10 and 20
res10_ca = None
res20_ca = None
for atom in atoms:
    if atom.name == 'CA':
        if atom.residue.id == 10:
            res10_ca = atom.index
        elif atom.residue.id == 20:
            res20_ca = atom.index

# Define collective variable (CV)
cv_force = CustomBondForce('r')
cv_force.addBond(res10_ca, res20_ca, [])
cv_force.setUsesPeriodicBoundaryConditions(True)

# Set up variables for metadynamics
cvs = [cv_force]
temperature = 300*kelvin
bias_factor = 10.0  # Commonly used value (5-20)
hill_height = 1.0*kilojoule_per_mole  # Height of the hills
hill_frequency = 100  # Add a hill every 100 steps

# Initialize metadynamics
meta = MetaDynamics(system, cvs, temperature, bias_factor, hill_height, hill_frequency)

# Prepare simulation
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Energy minimization
simulation.minimizeEnergy()

# Add reporters
simulation.reporters.append(DCDReporter('metadynamics_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('metadynamics_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# Custom reporter to record collective variable values
class CVReporter(object):
    def __init__(self, file, reportInterval):
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        self._file.write('Step,CV_Value,Bias_Energy\n')
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, True)
    
    def report(self, simulation, state):
        cv_value = meta.getCollectiveVariableValues(simulation.context)[0]
        bias_energy = meta.getBiasEnergy(simulation.context)
        self._file.write(f'{simulation.currentStep},{cv_value},{bias_energy}\n')

simulation.reporters.append(CVReporter('cv_values.csv', 100))

# Run simulation
simulation.step(500000)  # Run 500,000 steps of simulation (1 ns)
```

#### SimulatedTempering

```python
st = SimulatedTempering(system, temperatures, scaleForces=True)
```
- `system` (System) - System to simulate
- `temperatures` (list) - List of temperatures (with units)
- `scaleForces` (bool) - Whether to scale forces

Simulated tempering periodically changes the temperature during a simulation to allow crossing energy barriers. This enables exploration of conformational spaces that may not be accessible in regular constant-temperature simulations.

**Usage Example**:

```python
from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np

# Load structure from PDB file
pdb = PDBFile('protein.pdb')

# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Create basic integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Set temperature range (e.g., 8 levels between 300K and 400K)
min_temp = 300.0
max_temp = 400.0
num_temps = 8
temperatures = [min_temp + (max_temp - min_temp) * i / (num_temps-1) for i in range(num_temps)]
temperatures = [t*kelvin for t in temperatures]

# Initialize simulated tempering
st = SimulatedTempering(system, temperatures)

# Initialize free energy weights (initially equal, or use pre-calculated values)
weights = np.zeros(num_temps)
st.setTemperatureWeights(weights)

# Prepare simulation
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Energy minimization
simulation.minimizeEnergy()

# Add reporters
simulation.reporters.append(DCDReporter('simtemp_trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter('simtemp_output.csv', 1000, 
    step=True, time=True, potentialEnergy=True, temperature=True))

# Custom reporter to record temperature states
class TemperatureStateReporter(object):
    def __init__(self, file, reportInterval):
        self._file = open(file, 'w')
        self._reportInterval = reportInterval
        self._file.write('Step,TemperatureIndex,Temperature,AcceptanceRate\n')
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, True)
    
    def report(self, simulation, state):
        temp_index = st.getCurrentTemperature(simulation.context)
        temp = temperatures[temp_index].value_in_unit(kelvin)
        acceptance = st.getAcceptanceRate(temp_index)
        self._file.write(f'{simulation.currentStep},{temp_index},{temp},{acceptance}\n')

simulation.reporters.append(TemperatureStateReporter('temperature_states.csv', 100))

# Function to adaptively update weights
def update_weights(simulation, st, interval=50000):
    if simulation.currentStep % interval == 0 and simulation.currentStep > 0:
        counts = st.getStateVisitCounts()
        if min(counts) > 0:  # Ensure all states have been visited
            # Update weights based on state visit frequencies on a log scale
            weights = [-np.log(count/sum(counts)) for count in counts]
            # Set minimum to 0
            weights = [w - min(weights) for w in weights]
            st.setTemperatureWeights(weights)
            print(f"Step {simulation.currentStep}: Updated weights: {weights}")

# Run simulation (with adaptive weight updates)
for i in range(20):  # Divide total 1,000,000 steps (2 ns) into 20 parts
    simulation.step(50000)  # Run 50,000 steps at a time
    update_weights(simulation, st)  # Update weights
```

## Code Examples

### Example of Explicit Solvent System Setup

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# Load structure from PDB file
pdb = PDBFile('protein.pdb')

# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Solvate protein with water
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, padding=1.0*nanometer, model='tip3p')

# Build system
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Add temperature and pressure control (NPT ensemble)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# Set up integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Prepare simulation
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Energy minimization
simulation.minimizeEnergy()

# Equilibration
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(10000)  # 20 ps of equilibration

# Production simulation
simulation.reporters.append(DCDReporter('trajectory.dcd', 1000))
simulation.reporters.append(StateDataReporter(
    'data.csv', 1000, step=True, time=True,
    potentialEnergy=True, temperature=True, density=True
))
simulation.step(500000)  # 1 ns of production simulation
```

### Custom Force Field Example

```python
from openmm import *
from openmm.app import *
from openmm.unit import *

# Example of custom torsion potential
custom_torsion = CustomTorsionForce('k*(1+cos(n*theta-gamma))')
custom_torsion.addPerTorsionParameter('k')
custom_torsion.addPerTorsionParameter('n')
custom_torsion.addPerTorsionParameter('gamma')

# Example of custom nonbonded force (modified Lennard-Jones)
custom_nb = CustomNonbondedForce(
    '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)'
)
custom_nb.addPerParticleParameter('sigma')
custom_nb.addPerParticleParameter('epsilon')
```

## Troubleshooting

### Common Problems and Solutions

1. **Simulation is unstable and diverges**
   - Perform energy minimization properly
   - Reduce the time step (0.5-1 fs)
   - Review constraint conditions (especially bonds involving hydrogen atoms)
   - Check for steric clashes in the initial structure

2. **Performance issues**
   - Verify that an appropriate platform (CUDA/OpenCL) is being used
   - Set appropriate cutoffs for nonbonded interactions
   - Configure platform-specific optimization options

3. **Platform selection issues**
   - Explicitly select a specific platform:
     ```python
     platform = Platform.getPlatformByName('CUDA')
     properties = {'CudaPrecision': 'mixed'}
     simulation = Simulation(topology, system, integrator, platform, properties)
     ```

4. **Unit-related errors**
   - In OpenMM, physical quantities always need units
   - Use the `openmm.unit` module for unit conversion and calculation

## Function Reference

### openmm.unit Module

This module handles physical quantities with units.

```python
from openmm.unit import *

# Length units
nanometer      # Nanometer (1e-9 m)
angstrom       # Angstrom (1e-10 m)
picometer      # Picometer (1e-12 m)
bohr           # Bohr radius (quantum mechanical unit of length)

# Time units
picosecond     # Picosecond (1e-12 s)
femtosecond    # Femtosecond (1e-15 s)
nanosecond     # Nanosecond (1e-9 s)

# Energy units
kilojoule_per_mole  # Kilojoule/mole (biochemical standard)
kilocalorie_per_mole  # Kilocalorie/mole
hartree        # Hartree (quantum chemical energy unit)
electronvolt   # Electron volt

# Temperature units
kelvin         # Kelvin

# Mass units
dalton         # Dalton (atomic mass unit)
atomic_mass_unit  # Atomic mass unit

# Charge units
elementary_charge  # Elementary charge

# Pressure units
bar            # Bar
atmosphere     # Atmosphere

# Unit conversion
length = 5.0 * nanometer  # Length of 5 nm
length_in_angstrom = length.value_in_unit(angstrom)  # Convert to Angstroms
```

### Modeller Class

The `Modeller` class provides tools for modifying molecular systems.

```python
modeller = Modeller(topology, positions)
```
- `topology` (Topology) - Topology of the molecular system
- `positions` (list) - List of particle positions (with units)

#### Main Methods

**addSolvent**
```python
Modeller.addSolvent(forcefield, model='tip3p', boxSize=None, boxVectors=None,
                   padding=None, neutralize=True, positiveIon='Na+',
                   negativeIon='Cl-', ionicStrength=0*molar)
```
- `forcefield` (ForceField) - Force field to use
- `model` (str) - Water model ('tip3p', 'tip4pew', 'tip5p', etc.)
- `boxSize` (Vec3, optional) - Box size (with units)
- `boxVectors` (tuple, optional) - Periodic boundary condition box vectors
- `padding` (quantified length, optional) - Water padding around the solute (with units)
- `neutralize` (bool) - Whether to neutralize the system
- `positiveIon`, `negativeIon` (str) - Types of ions to use
- `ionicStrength` (quantified concentration) - Ionic strength (with units)

**addHydrogens**
```python
Modeller.addHydrogens(forcefield=None, pH=7.0, variants=None)
```
- `forcefield` (ForceField, optional) - Force field to use
- `pH` (float) - pH value
- `variants` (list, optional) - List of residue variants

**delete**
```python
Modeller.delete(atoms)
```
- `atoms` (list) - List of atoms to delete

### ForceField Class

The `ForceField` class handles the definition and application of molecular force fields.

```python
forcefield = ForceField(*files)
```
- `files` (str) - List of force field definition files

#### Main Methods

**createSystem**
```python
ForceField.createSystem(topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer,
                       constraints=None, rigidWater=True, removeCMMotion=True,
                       hydrogenMass=None, residueTemplates=dict(), verbose=False,
                       ewaldErrorTolerance=0.0005)
```
- `topology` (Topology) - Topology of the molecular system
- `nonbondedMethod` (int) - Method for calculating nonbonded interactions
- `nonbondedCutoff` (quantified length) - Cutoff distance for nonbonded interactions (with units)
- `constraints` (int) - Constraint level (`None`, `HBonds`, `AllBonds`, `HAngles`)
- `rigidWater` (bool) - Whether to treat water molecules as rigid
- `removeCMMotion` (bool) - Whether to remove center of mass motion
- `hydrogenMass` (quantified mass, optional) - Mass of hydrogen atoms (with units)
- `residueTemplates` (dict) - Mapping of residue templates
- `verbose` (bool) - Whether to display detailed output
- `ewaldErrorTolerance` (float) - Error tolerance for Ewald method

**getUnmatchedResidues**
```python
ForceField.getUnmatchedResidues(topology)
```
- `topology` (Topology) - Topology of the molecular system
- Returns: List of unmatched residues

**getMatchingTemplates**
```python
ForceField.getMatchingTemplates(topology, ignoreExternalBonds=False)
```
- `topology` (Topology) - Topology of the molecular system
- `ignoreExternalBonds` (bool) - Whether to ignore external bonds
- Returns: List of matching templates

### MonteCarloBarostat

```python
barostat = MonteCarloBarostat(pressure, temperature, frequency=25)
```
- `pressure` (quantified pressure) - Target pressure (with units)
- `temperature` (quantified temperature) - Simulation temperature (with units)
- `frequency` (int) - Frequency of Monte Carlo trials (number of steps)

### RMSDForce

```python
force = RMSDForce(referencePositions, particles=None)
```
- `referencePositions` (list) - List of reference positions (with units)
- `particles` (list, optional) - List of particle indices to include

Detailed documentation for OpenMM can be found on the [official website](http://docs.openmm.org/).