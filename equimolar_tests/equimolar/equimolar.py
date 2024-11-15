from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm
from ase.visualize import view
from multiprocessing import Pool
import os

n=200
N=1925
n_metals = len(metals)
composition = np.ones(n_metals)/n_metals

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)




def dissolve(i):
    np.random.seed(i)
    dissolver.make_particle(composition,n_atoms=N)
    
    atoms, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

    diss_comp = dissolver.get_composition(dissolved_atoms,precision=6)[1]

    
    N_final_111,final_111_comp = dissolver.get_composition_cn(atoms,9,precision=6)
    
    s = N_final_111/N_initial_111

    return np.concatenate((final_111_comp, diss_comp, [s]))

if __name__ == '__main__':
    with Pool(32) as pool:
        results = list(tqdm(pool.imap(dissolve,range(n)),total=n, desc= 'Processing'))

    results = np.array(results)
    final_111_compositions = results[:,:n_metals]
    dissolution_compositions = results[:,n_metals:-1]
    dissolution_factors = results[:,-1]


    data = np.hstack((final_111_compositions,dissolution_compositions,dissolution_factors.reshape(-1,1)))

    np.savetxt('equimolar/equimolar4nm_muM_results.csv',data,delimiter=',',header='111 composition, dissolution composition, dissoluton factor',fmt='%1.6f')
        




