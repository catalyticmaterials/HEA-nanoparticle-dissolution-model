from utilities.compositionspace_functions import get_molar_fractions
import numpy as np
from utilities.stability.particle_dissolver import Dissolver
from utilities.stability import metals,load_regressor
from multiprocessing import Pool
from tqdm import tqdm

mf = get_molar_fractions(0.1,5)

n=5
N=1289

model_path = '../regression_model/IrPdPtRhRu_multilinear.regressor'
regressor = load_regressor(model_path)

dissolver = Dissolver(metals)
dissolver.set_regressor(regressor)
dissolver.set_feature_generator('LinearNeighborID')
N_initial_111 = len(dissolver.get_111_atoms(dissolver.dummy_particle(N)))


def dissolve(comp):
    np.random.seed(42)
    fds=[]
    for _ in range(n):
        dissolver.make_particle(comp,N)
        atoms,_ = dissolver.dissolve_atoms(1.0,relax_cn=True)
        N_final_111 = len(dissolver.get_111_atoms(atoms))
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
    header='Ir,Pd,Pt,Rh,Ru,fd,fd_se'
    np.savetxt('grid_search_1V.csv',data,delimiter=',',fmt=fmt,header=header)