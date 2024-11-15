from utilities.particle_dissolver import Dissolver
from utilities import metals
from utilities.colors import metal_colors
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from ase.io import Trajectory

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)


def layer_profile(layer_compositions,natoms,ax):

    n_layers = layer_compositions.shape[0]

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
    

    return ax


    
def coordination_profile(cn_comps,ns,ax):


    x = np.arange(3,11)
    b=np.zeros(8)
    for i,metal in enumerate(metals):
        ax.bar(x,cn_comps[i],bottom=b,color=metal_colors[metal],width=0.9)
        b+=cn_comps[i]
    

    ax.set_yticks([0,0.2,0.4,0.6,0.8,1],[0,20,40,60,80,100])
    ax.set_xlabel('CN')
    ax.set_xticks(x,[i for i in range(3,10)] + ['>9'])
    ax.set_xlim(x[0]-0.55,x[-1]+0.55)
    ax.set_ylim(0,1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    for i,n in enumerate(ns):
        ax.text(x[i],1.0,f'{int(n)}',ha='center',va='bottom')

    return ax




final_layer_compositions = []
final_n_atoms_layer = []
cn_comps_all = np.zeros((8,8))
n_atoms_all = np.zeros(8)
for i in range(10):
    final_particle = Trajectory(f'trajectories/traj_{i}.traj')[-1]
    final_layer_compositions_, final_n_atoms_layer_ = dissolver.composition_by_layer(final_particle,n_layers=4)
    final_layer_compositions.append(final_layer_compositions_)
    final_n_atoms_layer.append(final_n_atoms_layer_)

    n_atoms = []
    cn_comps = []
    for cn in range(3,10):
        N_cn,comp_cn = dissolver.get_composition_cn(final_particle,cn,3)
        n_atoms.append(N_cn)
        cn_comps.append(comp_cn)

    symbols = np.array(final_particle.get_chemical_symbols())
    N_bulk,comp_bulk = dissolver.get_composition(symbols[dissolver.get_coordination_numbers(final_particle)>9])
    n_atoms.append(N_bulk)
    cn_comps.append(comp_bulk)

    cn_comps_all += np.array(cn_comps).T
    n_atoms_all += np.array(n_atoms)



cn_comps_all/=np.sum(cn_comps_all,axis=0)

final_layer_compositions = np.mean(final_layer_compositions,axis=0)
final_n_atoms_layer = np.sum(final_n_atoms_layer,axis=0)

fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(6,3),sharey=True,gridspec_kw={'width_ratios': [3, 5]})
layer_profile(final_layer_compositions,final_n_atoms_layer,ax1)

coordination_profile(cn_comps_all,n_atoms_all,ax2)

fig.legend(loc='outside upper center',ncols=8, mode='expand',shadow=False,fancybox=False)
plt.tight_layout()
plt.subplots_adjust(top=0.825)
plt.savefig('composition_profiles.png',dpi=600,bbox_inches='tight')

