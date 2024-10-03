from utilities.compositionspace_functions import get_random_molar_fractions
from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool


np.random.seed(2)
N=100   # composititons
n=30    # iterations

mfs = get_random_molar_fractions(8,N)

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)


def get_sd(mf):
    np.random.seed(42)
    Sds = []
    for i in range(n):
        dissolver.make_particle(mf)
        N_initial = len(dissolver.particle)

        atoms,diss = dissolver.dissolve_atoms(0.8,relax_cn=True)

        N_final = len(atoms)

        Sds.append(N_final/N_initial)

    Sd_data = [np.mean(Sds),np.std(Sds,ddof=1)]
    return np.append(mf,Sd_data)


if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(get_sd,mfs),total=N, desc= 'Processing'))


    header = ','.join(metals) + ',Sd_mean,Sd_std'
    np.savetxt('Sd_variance.csv',results,header=header,fmt='%1.6f',delimiter=',')