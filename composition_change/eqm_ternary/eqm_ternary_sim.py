from utilities.particle_dissolver import Dissolver
from utilities import metals
import numpy as np
from iteround import saferound
from tqdm import tqdm


dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

composition = saferound(np.array([1,1,1])/3,places=6)




for ternary_metals in (['Pd','Pt','Ru'],['Au','Cu','Pt'],['Au','Cu','Pd']):

    system = ''.join(ternary_metals)
    ternary_mask = np.array([metal in ternary_metals for metal in metals])

    full_composition = np.zeros(8)
    full_composition[ternary_mask] = composition


    initials = []
    finals = []
    diss = []
    Sd = []
    for i in tqdm(range(50),mininterval=10):
        np.random.seed(i)
        initial_particle = dissolver.make_particle(full_composition,n_atoms=1925,return_particle=True)
        N_initial_111,initial_111_comp = dissolver.get_composition_cn(initial_particle,9,4)

        final_particle, dissolved_atoms,traj_list = dissolver.dissolve_atoms_single(0.8,relax_func=dissolver.relax_particle_single_cn,return_trajectory=True)
        N_final_111,final_111_comp = dissolver.get_composition_cn(final_particle,9,4)

        traj_111_comp = []
        for atoms in traj_list:
            traj_111_comp.append(dissolver.get_composition_cn(atoms,9,4)[1])

        np.savetxt(f'eqm_ternary/{system}/trajectories/traj_111_comp_{i}.csv',np.array(traj_111_comp)[:,ternary_mask],delimiter=',',header=','.join(ternary_metals),fmt='%1.4f')


        initials.append(np.array(initial_111_comp)[ternary_mask])
        finals.append(np.array(final_111_comp)[ternary_mask])
        diss.append(np.array(dissolver.get_composition(dissolved_atoms)[1])[ternary_mask])
        Sd.append(N_final_111/N_initial_111)
        
    Sd = np.array(Sd).reshape(-1,1)
    header = '_i,'.join(ternary_metals) + '_i,' + '_f,'.join(ternary_metals) + '_f,' + '_d,'.join(ternary_metals) + '_d,Sd'
    np.savetxt(f'eqm_ternary/{system}/initial_and_final_111_comps.csv',np.hstack((initials,finals,diss,Sd)),delimiter=',',header=header,fmt='%1.4f')

