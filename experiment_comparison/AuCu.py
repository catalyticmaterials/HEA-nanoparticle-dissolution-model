from utilities.particle_dissolver import Dissolver
from ase.visualize import view
import numpy as np

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
np.random.seed(42)
dissolver.make_particle([0.,0.2,0.8,0.,0.,0.,0.,0.],n_atoms=4000)
view(dissolver.particle)
atoms,dissolved = dissolver.dissolve_atoms(0.8,relax_cn=True)

print(len(atoms))
view(atoms)
print(np.unique(dissolved,return_counts=True))