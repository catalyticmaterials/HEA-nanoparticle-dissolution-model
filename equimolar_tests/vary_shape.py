from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

n_metals = len(metals)
n=200
E_list = [1.0,1.1,1.2,1.3,1.4,1.5] # 100 energy relative to 111
N = 3355 # approx 5 nm
composition = np.ones(n_metals)/n_metals



#  Set up dissolver
dissolver = Dissolver(metals,regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

N_dict = {}
N_111_dict = {}
for E_100 in E_list:
    dummy = dissolver.dummy_particle(n_atoms=N,energies=[1,E_100])
    N_dict[str(E_100)] = len(dummy)
    N_111_dict[str(E_100)] = dissolver.get_composition_cn(dummy,9)[0]



def dissolve(E_100):
    # dummy = Dissolver.dummy_particle(n_atoms=N,energies=[1,E_100])
    # N_initial = len(dummy)
    # N_initial_111 = len(Dissolver.get_111_atoms(dummy))
    N_initial_111 = N_111_dict[str(E_100)]
    dissolver.make_particle(composition,n_atoms=N,energies=[1,E_100])
    atoms, dissolved_atoms = dissolver.dissolve_atoms(0.8,relax_cn=True)

    N_final_111, final_111_composition = dissolver.get_composition_cn(atoms,9)
    N_diss,dissolution_composition = dissolver.get_composition(dissolved_atoms,precision=6)


    dissolution_factor = N_final_111/N_initial_111
    return np.concatenate([final_111_composition,dissolution_composition,[dissolution_factor]])

    
if __name__ == '__main__':

    datalist = []
    for E_100 in E_list:
        Es = np.ones(n)*E_100
        np.random.seed(42)
        N_final=[]
        dissolution_compositions = []
        final_111_compositions = []
        dissolution_factors = []
        with Pool(4) as pool:
            results = list(tqdm(pool.imap(dissolve,Es),total=n, desc= 'Processing'))

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

        datalist.append(np.concatenate(([N_dict[str(E_100)],E_100],mean_111,std_111,mean_diss,std_diss,[mean_df,std_df])))

    header = 'N,E_100,mean_111,std_111,mean_diss,std_diss,mean_df,std_df'
    np.savetxt('shape_dependence.csv',datalist,delimiter=',',fmt='%1.6f',header=header)

