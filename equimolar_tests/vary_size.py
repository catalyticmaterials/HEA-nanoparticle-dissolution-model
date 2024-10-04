from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool



n=200
N_list = [201,405,711,1289,1925,3355] # approximate diameters (nm): 1.8, 2.3, 2.8, 3.5, 4.0, 5.0:
n_metals = len(metals)
composition = np.ones(n_metals)/n_metals

#  Set up dissolver
dissolver = Dissolver(metals,regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)


def dissolve(args):
    N,seed = args
    np.random.seed(seed)
    dissolver.make_particle(composition,n_atoms=N)
    N_initial_111,_ = dissolver.get_composition_cn(dissolver.particle,9)
    atoms, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

    N_final_111, final_111_composition = dissolver.get_composition_cn(atoms,9)
    N_diss,dissolution_composition = dissolver.get_composition(dissolved_atoms,precision=6)


    dissolution_factor = N_final_111/N_initial_111

    return np.concatenate([final_111_composition,dissolution_composition,[dissolution_factor]])


if __name__ == '__main__':

    datalist = []
    for N in N_list:
        args = [[N,i] for i in range(n)]
        N_final=[]
        dissolution_compositions = []
        final_111_compositions = []
        dissolution_factors = []
        with Pool(4) as pool:
            results = list(tqdm(pool.imap(dissolve,args),total=n, desc= 'Processing'))

        results = np.array(results)
        final_111_compositions = results[:,:n_metals]
        dissolution_compositions = results[:,n_metals:-1]
        dissolution_factors = results[:,-1]

        
        mean_111 = np.mean(final_111_compositions,axis=0)
        std_111 = np.std(final_111_compositions,axis=0,ddof=1)

        mean_diss = np.mean(dissolution_compositions,axis=0)
        std_diss = np.std(dissolution_compositions,axis=0,ddof=1)

        mean_df = np.mean(dissolution_factors)
        std_df = np.std(dissolution_factors,ddof=1)

        datalist.append(np.concatenate(([N],mean_111,std_111,mean_diss,std_diss,[mean_df,std_df])))

    header = 'N,mean_111,std_111,mean_diss,std_diss,mean_df,std_df'
    np.savetxt('size_dependence.csv',datalist,delimiter=',',fmt='%1.6f',header=header)

