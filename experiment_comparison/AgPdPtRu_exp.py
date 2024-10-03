from utilities.particle_dissolver import Dissolver
import numpy as np
from utilities import metals as all_metals
from multiprocessing import Pool
from tqdm import tqdm
import pandas as pd
from iteround import saferound




dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1)
N=6000
# n=5



def dissolve(full_comp):
    # for i in range(n):
    
    dissolver.make_particle(full_comp,N)
    atoms = dissolver.dissolve_atoms(0.85,relax_cn=True)[0]

    final_111_comp = dissolver.get_composition_cn(atoms,9,6)[1]

    return np.array(final_111_comp)
    

if __name__ == '__main__':

    for metals in (['Ag','Pt','Ru'],['Ag','Pd','Ru']):
        np.random.seed(42)
        system= ''.join(metals)
        
        metal_mask = np.array([metal in metals for metal in all_metals])

        df = pd.read_csv(f'{system}_Edx.csv',sep=',')
        comp = df[metals].to_numpy()
        
        neg_mask = comp<0
        comp[neg_mask] = 0.0
        comp = comp/np.sum(comp,axis=1).reshape(-1,1)
        comp = np.array([saferound(x,6) for x in comp])
        full_comp = np.zeros((len(comp),8))
        full_comp[:,metal_mask] = comp

        with Pool(4) as pool:
            results = list(tqdm(pool.imap(dissolve,full_comp),total=len(full_comp), desc= 'Processing'))


        data = np.hstack((comp,np.array(results)[:,metal_mask]))
        
        header='_i,'.join(metals) + '_i,' + '_f,'.join(metals) + '_f'
        np.savetxt(f'{system}_surface_comp_change.csv',data,delimiter=',',fmt='%1.6f',header=header)
            