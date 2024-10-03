import numpy as np
from utilities.stability.particle_dissolver import Dissolver
from utilities.stability import metals,load_regressor
from multiprocessing import Pool
from tqdm import tqdm

n=20
N=1289

model_path = '../regression_model/IrPdPtRhRu_multilinear.regressor'
regressor = load_regressor(model_path)

dissolver = Dissolver(metals)
dissolver.set_regressor(regressor)
dissolver.set_feature_generator('LinearNeighborID')
N_initial_111 = len(dissolver.get_111_atoms(dissolver.dummy_particle(N)))


def dissolve(U):
    np.random.seed(42)
    fds=[]
    for _ in range(n):
        dissolver.make_particle([0.2]*5,N)
        atoms,_ = dissolver.dissolve_atoms(U,relax_cn=False)
        N_final_111 = len(dissolver.get_111_atoms(atoms))
        fd = N_final_111/N_initial_111
        fds.append(fd)

    fd_mean = np.mean(fds)
    fd_se = np.std(fds,ddof=1)/np.sqrt(n)
    return fd_mean,fd_se


potentials=np.arange(0.0,1.26,0.01)
if __name__ == '__main__':
    with Pool(12) as pool:
        results = list(tqdm(pool.imap(dissolve,potentials),total=len(potentials), desc= 'Processing'))


    data = np.hstack((potentials.reshape(-1,1),results))
    fmt = ['%1.2f','%1.6f','%1.6f']
    header='U,fd,fd_se'
    np.savetxt('equimolar_potential_sweep_norelax.csv',data,delimiter=',',fmt=fmt,header=header)