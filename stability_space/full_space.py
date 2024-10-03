from utilities.compositionspace_functions import get_random_molar_fractions
from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm

np.random.seed(42)


mfs = get_random_molar_fractions(8,200)

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

data = []
for mf in tqdm(mfs):

    dissolver.make_particle(mf)
    N_initial = len(dissolver.particle)

    atoms,diss = dissolver.dissolve_atoms(0.8,relax_cn=True)

    N_final = len(atoms)

    Sd = N_final/N_initial

    data.append(np.append(mf,Sd))

header = ','.join(metals) + ',Sd'
np.savetxt('full_space_sampling.csv',data,header=header,fmt='%1.6f',delimiter=',')