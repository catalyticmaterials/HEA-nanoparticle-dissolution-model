from utilities.particle_dissolver import Dissolver
import numpy as np

N=1925

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
N_initial_111,_ = dissolver.get_composition_cn(dissolver.dummy_particle(N),9)


def dissolve(comp):

    data = []
    for i in range(10):

        np.random.seed(i)
        dissolver.make_particle(comp,n_atoms=N)
        
        atoms, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

        diss_comp = dissolver.get_composition(dissolved_atoms,precision=6)[1]

        
        N_final_111,final_111_comp = dissolver.get_composition_cn(atoms,9,precision=6)
        
        s = N_final_111/N_initial_111

        data.append(np.concatenate((final_111_comp, diss_comp, [s])))

    return np.mean(data,axis=0)