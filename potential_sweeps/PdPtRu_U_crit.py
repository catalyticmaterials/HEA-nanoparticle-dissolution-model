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
n = 5

ternary_metals = ['Pd','Pt','Ru']

metal_mask = np.array([metal in ternary_metals for metal in metals])


def dissolve_(comp):

    full_comp = np.zeros(8)
    full_comp[metal_mask] = comp

    np.random.seed(42)
    particle = dissolver.make_particle(full_comp,N,return_particle=True)

    # First take big potential steps
    for U in potentials:
        
        atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
        N_final_111,_ = dissolver.get_composition_cn(atoms,9)
   
        if N_final_111==0.0:
            # Go back to previous step and break
            U_min = U-U_step
            break

    U = U_min + U_fine_step
    # Increment U in smaller steps
    while True:
        atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
        N_final_111,_ = dissolver.get_composition_cn(atoms,9)

        if N_final_111==0:
            return U
        else:
            U+=U_fine_step


def dissolve(comp):
    

    full_comp = np.zeros(8)
    full_comp[metal_mask] = comp
    
    Uc = []
    for i in range(n):
        a = 0.0
        b = 2.0
        c = (b-a)/2

        U=a+c

        np.random.seed(i)
        particle = dissolver.make_particle(full_comp,N,return_particle=True)

        # First take big potential steps
        while (b-a)>0.005:
            
            
            atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
            N_final_111,_ = dissolver.get_composition_cn(atoms,9)
    
            if N_final_111==0.0:
                # Go back to previous step and break
                b=U
            else:
                a=U
            
            
            c=(b-a)/2
            U=a+c
        
        Uc.append(b)



        
    return Uc




mf = get_molar_fractions(0.05,3)

if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,mf),total=len(mf), desc= 'Processing'))


    data = np.hstack((mf,np.array(results)))
    fmt = ['%1.2f']*3+['%1.6f']*n
    header='Pd,Pt,Ru,U_crit'
    np.savetxt('PdPtRu_potential_sweep.csv',data,delimiter=',',fmt=fmt,header=header)