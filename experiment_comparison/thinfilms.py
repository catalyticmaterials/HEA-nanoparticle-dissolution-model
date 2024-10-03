from ase.build import fcc111
import numpy as np
from ase.visualize import view
from utilities.particle_dissolver import Dissolver

np.random.seed(42)


dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1)


applied_potential = 0.85

composition = np.array([0.3,0.0,0.0,0.0,0.0,0.3,0.0,0.4])

n=33
l=10
atoms = fcc111('Pt',size=(n,n,l),vacuum=10)
symbols = np.random.choice(dissolver.metals,size=len(atoms),p=composition)
# Set symbols
atoms.set_chemical_symbols(symbols)


# view(atoms)


# Simulation
while len(atoms)>(n*n):
    ids = np.arange(n*n,len(atoms))
    # if return_trajectory:
    #     trajectory_list.append(atoms.copy())

    # Get coordination numbers of all atoms
    coordination_numbers,nl = dissolver.get_coordination_numbers(atoms,return_neighborlist=True)
    coordination_numbers = coordination_numbers

    # Symbol and IDs of atoms
    symbols = np.array(atoms.get_chemical_symbols())

    
    # Mask of atoms on surface to calculate
    surface_mask = (coordination_numbers[ids]<=dissolver.max_cn)*(coordination_numbers[ids]>=dissolver.min_cn)
    if np.sum(surface_mask)>0:
        surface_ids = ids[surface_mask]
        surface_symbols = np.array(symbols)[ids][surface_mask]

        # get featuers of surface atoms
        features = dissolver.feature_generator(nl,surface_ids,symbols,coordination_numbers)

        # Get dissolution potentials
        dissolution_potentials = dissolver.get_dissolution_potential(surface_symbols,features)

        # Dissolve atoms where U_diss<U
        remove_ids = surface_ids[dissolution_potentials<applied_potential]
    else:
        # Empty array if there are no atoms on surface with specified CN range
        remove_ids = np.empty(0,dtype=int)

    # Dissolve atoms with lower CN than included in model
    if np.any(coordination_numbers<dissolver.min_cn):
        remove_ids = np.append(remove_ids,ids[coordination_numbers[ids]<dissolver.min_cn])

    if len(remove_ids)==0:
        # Break simulation if there are no atoms to dissolve
        break
    else:
        # List the atoms to dissolve's position as free
        # free_positions= np.vstack((free_positions,np.atleast_2d(atoms.positions[remove_ids])))
        # Save symbols of atoms to dissolve
        # removed_atoms = np.append(removed_atoms,symbols[remove_ids])
        # Dissolve atoms by deleting them from ASE atoms object
        del atoms[remove_ids]


ids = np.arange(n*n,len(atoms))
cn = dissolver.get_coordination_numbers(atoms)[ids]
symbols = np.array(atoms.get_chemical_symbols())[ids]

print(dissolver.get_composition(symbols[cn<12]))
print(dissolver.get_composition(symbols[cn==9]))

view(atoms)