import numpy as np
from utilities.particle_dissolver import Dissolver
from multiprocessing import Pool
from tqdm import tqdm

n=30
N=1289

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)

composition = [0.2,0.2,0.2,0.0,0.2,0.2,0.0,0.0]

def dissolve(U):
    np.random.seed(42)
    fds=[]
    for _ in range(n):
        dissolver.make_particle(composition,N)
        atoms,_ = dissolver.dissolve_atoms(U,relax_cn=True)
        N_final_111,_ = dissolver.get_composition_cn(atoms,9)
        fd = N_final_111/N_initial_111
        fds.append(fd)

    fd_mean = np.mean(fds)
    std = np.std(fds,ddof=1)
    fd_se = std/np.sqrt(n)
    return fd_mean,fd_se, std


potentials=np.arange(-0.5,1.51,0.05)
if __name__ == '__main__':
    with Pool(4) as pool:
        results = list(tqdm(pool.imap(dissolve,potentials),total=len(potentials), desc= 'Processing'))


    data = np.hstack((potentials.reshape(-1,1),results))
    fmt = ['%1.2f','%1.6f','%1.6f','%1.6f']
    header='U,fd,fd_se,std'
    np.savetxt('equimolar_GPGM_potential_sweep.csv',data,delimiter=',',fmt=fmt,header=header)