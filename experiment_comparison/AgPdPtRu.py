from utilities.particle_dissolver import Dissolver
from utilities.compositionspace_functions import make_ternary_plot, get_molar_fractions
import numpy as np
from utilities import metals as all_metals
from multiprocessing import Pool
from tqdm import tqdm

mf = get_molar_fractions(0.1,4)
AgPdRu_mask = mf[:,2]==0
AgPtRu_mask = mf[:,1]==0
ternary_mask =  AgPdRu_mask + AgPtRu_mask
ternary_comps = mf[ternary_mask]

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
N=6000
# n=5
N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)

metals = ['Ag','Pd','Pt','Ru']

metal_mask = np.array([metal in metals for metal in all_metals])

def dissolve(comp):
    full_comp = np.zeros(8)
    full_comp[metal_mask] = comp
    # for i in range(n):
    np.random.seed(42)
    dissolver.make_particle(full_comp,N)
    atoms,diss = dissolver.dissolve_atoms(0.85,relax_cn=True)

    N_111, final_111_comp = dissolver.get_composition_cn(atoms,9)
    layer_comp,layer_natoms = dissolver.composition_by_layer(atoms,n_layers=4)
    if np.sum(layer_natoms[:3])==0:
        final_3layer_comp = np.zeros(8)
    else:
        final_3layer_comp = np.average(layer_comp[:3],axis=0,weights=layer_natoms[:3])
    Sd = [N_111/N_initial_111]

    return np.concatenate((comp,np.array(final_111_comp)[metal_mask],final_3layer_comp[metal_mask],Sd))
    

if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,mf),total=len(mf), desc= 'Processing'))


    data = np.array(results)
    fmt = ['%1.2f']*12+['%1.4f']
    header='Ag_b,Pd_b,Pt_b,Ru_b,Ag_f,Pd_f,Pt_f,Ru_f,Ag_f3,Pd_f3,Pt_f3,Ru_f3,S_d'
    np.savetxt('AgPdPtRu.csv',data,delimiter=',',fmt=fmt,header=header)
        