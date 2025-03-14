from utilities.particle_dissolver import Dissolver
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

N=3355
n=50

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
composition = np.ones(8)/8
dummy = dissolver.dummy_particle(N)
N_initial,_ = dissolver.get_composition_cn(dummy,9)




def run_simulation(args):

    dissolve_func,relax_func,particle = args  
    
    np.random.seed(42)
    atoms,dissolved = dissolve_func(0.8,atoms=particle,relax_func=relax_func)

    (surf_comp,bulk_comp),(n_atoms) = dissolver.composition_by_layer(atoms,n_layers=2,precision=3)


    N_111,composition_111 = dissolver.get_composition_cn(atoms,9,3)

    dissolved_comp = dissolver.get_composition(dissolved,3)[1]

    Sd = N_111/N_initial

    CNs = dissolver.get_coordination_numbers(atoms)
    cn_dist = [np.sum(CNs==cn) for cn in range(3,13)]

    return np.concatenate((surf_comp,composition_111,bulk_comp,dissolved_comp,n_atoms,cn_dist,[Sd,len(atoms)]))




if __name__ == '__main__':

    particles = []
    for i in range(n):
        np.random.seed(i)
        particles.append(dissolver.make_particle(composition,N,return_particle=True))
    
    for name in ('single_relax','batch_relax','batch_none'):

        if name.startswith('batch'):
            diss_func = dissolver.dissolve_atoms_batch
            if name.endswith('none'):
                relax_func=None
            else:
                relax_func=dissolver.relax_particle_batch_cn
        else:
            diss_func = dissolver.dissolve_atoms_single
            relax_func = dissolver.relax_particle_single_cn

        
        args = [[diss_func,relax_func,particle] for particle in particles]

        with Pool(4) as pool:
            results = list(tqdm(pool.imap(run_simulation,args),total=n, desc= 'Processing'))

        

        np.savetxt(f'model_variations/{name}.csv',results,delimiter=',',fmt='%1.3f')
        