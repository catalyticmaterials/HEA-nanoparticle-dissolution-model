from utilities.particle_dissolver import Dissolver
from ase.visualize import view
import numpy as np

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1)
np.random.seed(42)
dissolver.make_particle([0.0,0.0,0.25,0.0,0.0,0.75,0.0,0.0])
view(dissolver.particle)
atoms,dissolved = dissolver.dissolve_atoms(0.8,relax_cn=True)


view(atoms)

