import numpy as np
from utilities.particle_dissolver import Dissolver
from utilities import metals
from tqdm import tqdm
from multiprocessing import Pool

data  = np.loadtxt('grid_sim/grids_full.csv',delimiter=',',usecols=(0,1,2,3,4,5,6,7,24))

N=1925
dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)


def add_isc(dataline):
    comp = dataline[:8]
    rest = dataline[8:]
    isc = []    
    for i in range(10):
        np.random.seed(i)
        initial_particle = dissolver.make_particle(comp,n_atoms=N,return_particle=True)
        N_initial_111,initial_111_comp = dissolver.get_composition_cn(initial_particle,9,precision=6)
        isc.append(initial_111_comp)
    isc_mean = np.mean(isc,axis=0)
    return np.concatenate((comp,isc_mean,rest))


if __name__ == '__main__':
    with Pool(32) as pool:
        results = list(tqdm(pool.imap(add_isc,data),total=len(data), desc= 'Processing',mininterval=1))
    header = ','.join(metals) + ',' + '_isc,'.join(metals) + '_isc,' + '_fsc,'.join(metals) + '_fsc,' + '_d,'.join(metals) + '_d.Sd'
    np.savetxt('grid_data.csv',results,delimiter=',',fmt='%1.6f',header=header)