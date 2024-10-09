from utilities.particle_dissolver import Dissolver
from utilities import metals
from utilities.colors import metal_colors
import numpy as np
import matplotlib.pyplot as plt
from ase.visualize import view
from matplotlib.ticker import MultipleLocator

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)


def layer_profile(layer_compositions,natoms,ax):


    
    n_layers = layer_compositions.shape[0]


    # fig,ax = plt.subplots(figsize=(4,2.5))

    x = range(1,n_layers+1)
    b=np.zeros(n_layers)
    for i,metal in enumerate(metals):
        ax.bar(x,layer_compositions.T[i],bottom=b,color=metal_colors[metal],width=0.9,label=metal)
        b+=layer_compositions.T[i]


    ax.set_yticks([0,0.2,0.4,0.6,0.8,1],[0,20,40,60,80,100])
    ax.set_ylabel('Composition (at. %)')
    ax.set_xlabel('Surface Layer')
    ax.set_xticks(np.arange(1,len(natoms)+1),[i for i in np.arange(1,len(natoms))] + ['Bulk'])
    ax.set_xlim(1-0.55,len(natoms)+0.55)
    ax.set_ylim(0,1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    for i,n in enumerate(natoms):
        ax.text(x[i],1.0,f'{n}',ha='center',va='bottom')
    
    # plt.subplots_adjust(right=0.8)
    # fig.legend(loc='outside center right',ncols=1)

    return ax


def coordination_profile(atoms,ax):
    # fig,ax = plt.subplots(figsize=(4,2.5))

    n_atoms = []
    cn_comps = []
    for cn in range(3,10):
        N_cn,comp_cn = dissolver.get_composition_cn(atoms,cn,3)
        n_atoms.append(N_cn)
        cn_comps.append(comp_cn)

    symbols = np.array(atoms.get_chemical_symbols())
    N_bulk,comp_bulk = dissolver.get_composition(symbols[dissolver.get_coordination_numbers(atoms)>9])
    n_atoms.append(N_bulk)
    cn_comps.append(comp_bulk)

    cn_comps = np.array(cn_comps).T


    x = np.arange(3,11)
    b=np.zeros(8)
    for i,metal in enumerate(metals):
        ax.bar(x,cn_comps[i],bottom=b,color=metal_colors[metal],width=0.9)
        b+=cn_comps[i]
    

    ax.set_yticks([0,0.2,0.4,0.6,0.8,1],[0,20,40,60,80,100])
    # ax.set_ylabel('Composition (at. %)')
    ax.set_xlabel('CN')
    ax.set_xticks(x,[i for i in range(3,10)] + ['>9'])
    ax.set_xlim(x[0]-0.55,x[-1]+0.55)
    ax.set_ylim(0,1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    for i,n in enumerate(n_atoms):
        ax.text(x[i],1.0,f'{n}',ha='center',va='bottom')
    
    # plt.subplots_adjust(right=0.8)
    # fig.legend(loc='outside center right',ncols=1)
    return ax
    
        

np.random.seed(2)

# 8 metallic
composition = np.ones(8)/8

dissolver.make_particle(composition,1925)
initial_layer_compositions, initial_n_atoms_layer = dissolver.composition_by_layer(dissolver.particle,n_layers=5)
print(initial_layer_compositions)

final_particle,diss = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,traj_file='layer_profiles/eqm_all_metals.traj')

final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(6,3),sharey=True)
layer_profile(final_layer_compositions,final_n_atoms_layer,ax1)

# plt.tight_layout()
# plt.savefig('layer_profiles/eqm_all_metals_final.png',dpi=600,bbox_inches='tight')
# plt.close()
print(dissolver.get_composition(diss))
coordination_profile(final_particle,ax2)
# plt.savefig('layer_profiles/eqm_all_metals_final_CN.png',dpi=600,bbox_inches='tight')
# plt.close()
fig.legend(loc='outside upper center',ncols=8, mode='expand',shadow=False,fancybox=False)
plt.tight_layout()
plt.subplots_adjust(top=0.825)
plt.savefig('layer_profiles/eqm_all_metals.png',dpi=600,bbox_inches='tight')
print(final_layer_compositions)


stop

np.random.seed(1)

# PGM
composition = np.array([0.0]*3 + [0.2]*5)

dissolver.make_particle(composition,1925)

final_particle,diss = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

layer_profile(final_layer_compositions,final_n_atoms_layer)

# plt.tight_layout()
plt.savefig('layer_profiles/eqm_pgm_final.png',dpi=600,bbox_inches='tight')
plt.close()

coordination_profile(final_particle)
plt.savefig('layer_profiles/eqm_pgm_final_CN.png',dpi=600,bbox_inches='tight')
plt.close()





# GPGM
composition = np.array([0.2,0.2,0.2,0.0,0.2,0.2,0.0,0.0])

dissolver.make_particle(composition,1925)

final_particle,diss = dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn)

final_layer_compositions, final_n_atoms_layer = dissolver.composition_by_layer(final_particle,n_layers=5)

layer_profile(final_layer_compositions,final_n_atoms_layer)

# plt.tight_layout()
plt.savefig('layer_profiles/eqm_gpgm_final.png',dpi=600,bbox_inches='tight')
plt.close()

coordination_profile(final_particle)
plt.savefig('layer_profiles/eqm_gpgm_final_CN.png',dpi=600,bbox_inches='tight')
plt.close()