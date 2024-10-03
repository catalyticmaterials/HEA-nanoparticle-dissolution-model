from utilities.compositionspace_functions import get_molar_fractions
import numpy as np
from utilities.particle_dissolver import Dissolver
from utilities import metals
from multiprocessing import Pool
from tqdm import tqdm

mf = get_molar_fractions(0.1,8)


n=5
N=1289



dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)


def dissolve(comp):
    np.random.seed(42)
    fds=[]
    for _ in range(n):
        dissolver.make_particle(comp,N)
        atoms,_ = dissolver.dissolve_atoms(0.8,relax_cn=True)
        N_final_111,_ = dissolver.get_composition_cn(atoms,9)
        fd = N_final_111/N_initial_111
        fds.append(fd)

    fd_mean = np.mean(fds)
    fd_se = np.std(fds,ddof=1)/np.sqrt(n)
    return fd_mean,fd_se



if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,mf),total=len(mf), desc= 'Processing'))


    data = np.hstack((mf,results))
    fmt = ['%1.2f']*5 + ['%1.6f','%1.6f']
    header='Ag,Au,Cu,Ir,Pd,Pt,Rh,Ru,fd,fd_se'
    np.savetxt('grid_search_08V.csv',data,delimiter=',',fmt=fmt,header=header)