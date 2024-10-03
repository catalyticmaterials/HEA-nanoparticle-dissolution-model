from utilities.compositionspace_functions import get_molar_fractions
from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm




mfs = get_molar_fractions(0.1,5)

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

for quinary_metals,name in zip((['Ag','Au','Cu','Pd','Pt'],['Ir','Pd','Pt','Rh','Ru']),('GPGM','PGM')):

    metal_mask = np.array([metal in quinary_metals for metal in metals])

    data = []
    for mf in tqdm(mfs):

        full_mf = np.zeros(8)
        full_mf[metal_mask] = mf
        np.random.seed(42)
        dissolver.make_particle(full_mf)
        N_initial = len(dissolver.particle)

        atoms,diss = dissolver.dissolve_atoms(0.8,relax_cn=True)

        N_final = len(atoms)

        Sd = N_final/N_initial

        data.append(np.append(mf,Sd))

    header = ','.join(quinary_metals) + ',Sd'
    np.savetxt(f'{name}_grid.csv',data,header=header,fmt='%1.6f',delimiter=',')