from utilities.particle_dissolver import Dissolver
from utilities import metals
from utilities.colors import metal_colors
import numpy as np
import matplotlib.pyplot as plt
from ase.visualize import view
from ase.db import connect
from ase.io import write



dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

binary_metals = ['Cu','Pt']

binary_composition = np.array([0.25,0.75])

composition = np.zeros(8)
metal_mask = np.array([metal in binary_metals for metal in metals])
composition[metal_mask] = binary_composition





def binary_layer_profile(layer_compositions,natoms,m1_initial,m1,m2):

    M1 = layer_compositions[:,metals.index(m1)]
    M2 = layer_compositions[:,metals.index(m2)]
    
    
    se = np.sqrt(M1*M2/natoms)

    fig,ax = plt.subplots(figsize=(4,3))

    x = range(1,len(M1)+1)

    ax.bar(x,M1,color=metal_colors[m1],width=0.9,label=m1)
    ax.bar(x,M2,bottom=M1,color=metal_colors[m2],width=0.9,label=m2)

    ax.errorbar(x,M1,se,c='k',fmt='none',capsize=3,linewidth=1)

    ax.axhline(m1_initial,color='k',ls=':')

    ax.set_yticks([0,0.2,0.4,0.6,0.8,1],[0,20,40,60,80,100])
    ax.set_ylabel('Composition (at. %)')
    ax.set_xlabel('Atomic Surface Layer')
    ax.set_xticks(np.arange(1,len(natoms)+1),[i for i in np.arange(1,len(natoms))] + ['Bulk'])
    ax.set_xlim(1-0.55,len(natoms)+0.55)
    ax.set_ylim(0,1)

    for i,n in enumerate(natoms):
        ax.text(x[i],1.0,f'{n}',ha='center',va='bottom')
    
    ax.legend(loc='lower right',ncols=2)

    return fig,ax



np.random.seed(2)
initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)

initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(initial_particle,n_layers=5)


binary_layer_profile(initial_layer_compositions,initial_n_atoms_layer,0.75,'Pt','Cu')
plt.savefig('layer_profiles/Pt75Cu25_initial.png',dpi=600,bbox_inches='tight')



final_particle, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/Pt75Cu25.traj')
final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)


binary_layer_profile(final_layer_compositions,final_n_atoms_layer,0.75,'Pt','Cu')
plt.savefig('layer_profiles/Pt75Cu25_final.png',dpi=600,bbox_inches='tight')




# PtCu3

binary_composition = np.array([0.75,0.25])

composition[metal_mask] = binary_composition

np.random.seed(2)
initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)

initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(initial_particle,n_layers=5)


binary_layer_profile(initial_layer_compositions,initial_n_atoms_layer,0.75,'Cu','Pt')
plt.savefig('layer_profiles/Pt25Cu75_initial.png',dpi=600,bbox_inches='tight')



final_particle, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/Pt25Cu75.traj')
final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)


binary_layer_profile(final_layer_compositions,final_n_atoms_layer,0.75,'Cu','Pt')
plt.savefig('layer_profiles/Pt25Cu75_final.png',dpi=600,bbox_inches='tight')





# PtRu

binary_metals = ['Pt','Ru']
metal_mask = np.array([metal in binary_metals for metal in metals])
binary_composition = np.array([0.25,0.75])
composition = np.zeros(8)
composition[metal_mask] = binary_composition

np.random.seed(2)
initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)

initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(initial_particle,n_layers=5)


binary_layer_profile(initial_layer_compositions,initial_n_atoms_layer,0.75,'Ru','Pt')
plt.savefig('layer_profiles/PtRu_initial.png',dpi=600,bbox_inches='tight')



final_particle, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/PtRu.traj')
final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

print(dissolver.get_composition(final_particle.get_chemical_symbols())[1])

binary_layer_profile(final_layer_compositions,final_n_atoms_layer,0.75,'Ru','Pt')
plt.savefig('layer_profiles/PtRu_final.png',dpi=600,bbox_inches='tight')






# PdRu
binary_metals = ['Pd','Ru']
metal_mask = np.array([metal in binary_metals for metal in metals])
binary_composition = np.array([0.75,0.25])
composition = np.zeros(8)
composition[metal_mask] = binary_composition

np.random.seed(2)
initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)

initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(initial_particle,n_layers=5)


binary_layer_profile(initial_layer_compositions,initial_n_atoms_layer,0.75,'Pd','Ru')
plt.savefig('layer_profiles/PdRu_initial.png',dpi=600,bbox_inches='tight')



final_particle, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/PdRu.traj')
final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

print(dissolver.get_composition(final_particle.get_chemical_symbols())[1])

binary_layer_profile(final_layer_compositions,final_n_atoms_layer,0.75,'Pd','Ru')
plt.savefig('layer_profiles/PdRu_final.png',dpi=600,bbox_inches='tight')





# AuPd
binary_metals = ['Au','Pd']
metal_mask = np.array([metal in binary_metals for metal in metals])
binary_composition = np.array([0.25,0.75])
composition = np.zeros(8)
composition[metal_mask] = binary_composition

np.random.seed(2)
initial_particle = dissolver.make_particle(composition,n_atoms=1925,return_particle=True)

initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(initial_particle,n_layers=5)


binary_layer_profile(initial_layer_compositions,initial_n_atoms_layer,0.75,'Pd','Au')
plt.savefig('layer_profiles/AuPd_initial.png',dpi=600,bbox_inches='tight')



final_particle, dissolved_atoms = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/AuPd.traj')
final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

print(dissolver.get_composition(final_particle.get_chemical_symbols())[1])

binary_layer_profile(final_layer_compositions,final_n_atoms_layer,0.75,'Pd','Au')
plt.savefig('layer_profiles/AuPd_final.png',dpi=600,bbox_inches='tight')