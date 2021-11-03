# Use the OpenFF Toolkit to generate a minimal chemical topology
from openff.toolkit.topology import Molecule
from openff.toolkit.utils.toolkits import ToolkitRegistry, AmberToolsToolkitWrapper

#molecule = Molecule.from_smiles("CCO")
#molecule = Molecule.from_file("styrene.sdf")
#molecule =  Molecule.from_smiles("CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\CC/C=C\C")
#molecule =  Molecule.from_smiles("c1ccccc1")

toolkit_registry = ToolkitRegistry()
toolkit_registry.register_toolkit(AmberToolsToolkitWrapper)

molname = "molecules/br20.sdf"
molecule = Molecule.from_file(molname)
#molecule = Molecule.from_smiles('CCCCCC')
molecule.generate_conformers()

#print(molecule.conformers)
print('calculating am1-bcc...')
#am1_bcc_charge = 
#toolkit_registry.call('compute_partial_charges_am1bcc', molecule) # return charges::numpy.array of shape (natoms) of type float
molecule.compute_partial_charges_am1bcc(use_conformers=molecule.conformers, toolkit_registry=toolkit_registry)
print('done!')
exit()

#molecule.generate_conformers(n_conformers=1)
topology = molecule.to_topology()

# Load OpenFF 2.0.0 "Sage"
from openff.toolkit.typing.engines.smirnoff import ForceField
sage = ForceField("openff-2.0.0.offxml")

# Create an Interchange object
from openff.interchange.components.interchange import Interchange
out = Interchange.from_smirnoff(force_field=sage, topology=topology)

# Define box vectors and assign atomic positions
import numpy as np
out.box = [40, 40, 40] * np.eye(3)
out.positions = molecule.conformers[0]

# Convert the Interchnage object to an OpenMM System
omm_sys = out.to_openmm(combine_nonbonded_forces=True)

# or write to GROMACS files
#out.to_gro("out.gro")
#out.to_top("out.top")

#print(topology.charge_model)
out.to_lammps( molname.replace(".sdf", ".data") )

#from openff.interchange.drivers.lammps import _write_lammps_input
#_write_lammps_input(out, "run.inp")