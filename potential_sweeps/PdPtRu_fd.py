import numpy as np
from utilities.particle_dissolver import Dissolver
from utilities import metals
from utilities.compositionspace_functions import get_molar_fractions
from multiprocessing import Pool
from tqdm import tqdm

N=1289

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)

U_step = 0.1
U_fine_step = 0.01
potentials=np.arange(0.3,2.0,U_step)

ternary_metals = ['Pd','Pt','Ru']

metal_mask = np.array([metal in ternary_metals for metal in metals])

potentials=np.arange(0.2,1.1,0.01)


def dissolve(comp):

    full_comp = np.zeros(8)
    full_comp[metal_mask] = comp

    np.random.seed(42)
    particle = dissolver.make_particle(full_comp,N,return_particle=True)

    fds = []
    # First take big potential steps
    for U in potentials:
        
        atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
        N_final_111,_ = dissolver.get_composition_cn(atoms,9)
        fds.append(N_final_111/N_initial_111)

    return fds
        
        




mf = get_molar_fractions(0.05,3)

if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,mf),total=len(mf), desc= 'Processing'))


    data = np.hstack((mf,np.array(results)))
    header='Pd,Pt,Ru,fd(U)'
    np.savetxt('PdPtRu_potential_sweep_fd.csv',data,delimiter=',',fmt='%1.4f',header=header)