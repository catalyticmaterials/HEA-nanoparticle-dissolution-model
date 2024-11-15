from utilities.particle_dissolver import Dissolver
import numpy as np
from utilities import metals
from tqdm import tqdm

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
n=10 # number of particles for each composition

N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(1925),9)


metal_mask = [metal in ['Au','Pd'] for metal in metals]
data = []
for Au in tqdm(np.arange(0,1.05,0.05)):

    Pd = 1-Au

    composition = np.array([0.0,Au,0.0,0.0,Pd,0.0,0.0,0.0])
    Sds = []
    final_comps = []
    for i in range(n):

        np.random.seed(i)
        
        dissolver.make_particle(composition,1925)

        atoms,_=dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

        N_final_111, final_111_comp = dissolver.get_composition_cn(atoms,9,4)

        Sds.append(N_final_111/N_initial_111)
        final_comps.append(final_111_comp)
    
    Sd_mean = np.mean(Sds)
    Sd_se = np.std(Sds,ddof=1)/np.sqrt(n)
    mean_finals_comp = np.mean(final_comps,axis=0)[metal_mask]
    data.append(np.concatenate(([Au,Pd],mean_finals_comp,[Sd_mean,Sd_se])))


np.savetxt('Au_protection/PdAu.csv',data,delimiter=',',header='Au_i,Pd_i,Au_f,Pd_f,Sd_mu,Sd_se',fmt='%1.4f')
