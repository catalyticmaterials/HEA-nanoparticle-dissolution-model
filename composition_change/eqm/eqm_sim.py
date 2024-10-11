from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from iteround import saferound
from tqdm import tqdm
from multiprocessing import Pool


dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

composition = np.ones(8)/8
n=10

def dissolve(i):

    np.random.seed(i)
    initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)
    N_initial_111,initial_111_comp = dissolver.get_composition_cn(initial_particle,9,4)

    final_particle, dissolved_atoms,traj_list = dissolver.dissolve_atoms_single(0.8,relax_func=dissolver.relax_particle_single_cn,return_trajectory=True)
    N_final_111,final_111_comp = dissolver.get_composition_cn(final_particle,9,4)

    traj_111_comp = []
    for atoms in traj_list:
        traj_111_comp.append(dissolver.get_composition_cn(atoms,9,4)[1])
    np.savetxt(f'eqm/trajectories/traj_111_comp_{i}.csv',np.array(traj_111_comp),delimiter=',',header=','.join(metals),fmt='%1.4f')

    diss_comp = dissolver.get_composition(dissolved_atoms,4)[1]

    Sd = [N_final_111/N_initial_111]

    return np.concatenate((initial_111_comp,final_111_comp,diss_comp,Sd))


if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,range(n)),total=n, desc= 'Processing'))        

    header = '_i,'.join(metals) + '_i,' + '_f,'.join(metals) + '_f,' + '_d,'.join(metals) + '_d,Sd'
    np.savetxt(f'eqm/initial_and_final_111_comps.csv',results,delimiter=',',header=header,fmt='%1.4f')

