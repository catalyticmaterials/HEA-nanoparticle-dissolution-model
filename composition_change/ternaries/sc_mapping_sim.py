from utilities import metals
from utilities.particle_dissolver import Dissolver
from utilities.compositionspace_functions import get_molar_fractions
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)




compositions = get_molar_fractions(0.05,3)

def dissolve(composition):
    sds = []
    iscs = []
    fscs = []
    for i in range(10):
        np.random.seed(i)
        initial_particle = dissolver.make_particle(composition,n_atoms=1289,return_particle=True)
        N_initial_111,initial_111_comp = dissolver.get_composition_cn(initial_particle,9,precision=3)
        
        final_particle, dissolved_atoms = dissolver.dissolve_atoms(0.8,relax_cn=True)
        N_final_111,final_111_comp = dissolver.get_composition_cn(final_particle,9,precision=3)
        
        sds.append(N_final_111/N_initial_111)
        iscs.append(initial_111_comp)
        fscs.append(final_111_comp)
    fscs_mean = np.mean(fscs,axis=0)
    sum_ = np.sum(fscs_mean)
    if sum_!=0.0 and sum_!=1.0:
        fscs_mean = fscs_mean/sum_
    return np.concatenate((composition,np.mean(iscs,axis=0),fscs_mean,[np.mean(sds)]))




if __name__ == '__main__':
    for ternary_metals in (['Pd','Pt','Ru'],['Au','Cu','Pt'],['Au','Cu','Pd']):

        system = ''.join(ternary_metals)
  

        ternary_mask = [metal in ternary_metals for metal in metals]

        full_compositions = np.zeros((compositions.shape[0],len(metals)))
        full_compositions[:,ternary_mask] = compositions
    
        with Pool(4) as pool:
            results = list(tqdm(pool.imap(dissolve,full_compositions),total=len(full_compositions), desc= 'Processing'))


        data_mask = ternary_mask*3 + [True]

        data = np.array(results)[:,data_mask]

        header = '_b,'.join(ternary_metals) + '_b,' + '_isc,'.join(ternary_metals) + '_isc,' + '_fsc,'.join(ternary_metals) + '_fsc,Sd'
        np.savetxt(f'ternaries/{system}_sc_maping.csv',data,header=header,fmt='%1.3f',delimiter=',')

 