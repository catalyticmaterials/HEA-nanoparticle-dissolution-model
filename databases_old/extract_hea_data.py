from ase.db import connect
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from utilities.stability import metals, U_standard
from utilities.colors import metal_colors



with connect('bulks.db') as db:
    E_bulk={row.metal:row.energy for row in db.select()}


# Number of different slab compositions
N = 50



# Function to handle 'no match' error in database
def util_func(db,slab_idx,metal):
    try: 
        return db.get(slab_idx=slab_idx,metal=metal,defect='replace').energy
    except KeyError:
        return db.get(slab_idx=slab_idx,metal=metal,defect='none').energy


def ax_violin(ax1,ax2,dE,U,cn,metal):
    parts1 = ax1.violinplot(dE,[cn],showextrema=False)
    parts2 = ax2.violinplot(U,[cn],showextrema=False)

    for pc in parts1['bodies']:
        pc.set_facecolor(metal_colors[metal])
        pc.set_edgecolor(metal_colors[metal])
        pc.set_alpha(1)

    for pc in parts2['bodies']:
        pc.set_facecolor(metal_colors[metal])
        pc.set_edgecolor(metal_colors[metal])
        pc.set_alpha(1)

    



def get_nearest_neighbors(atoms,idx):
    cutoff = natural_cutoffs(atoms,mult=0.9)


    nl = NeighborList(cutoff,self_interaction=False, bothways=True)

    # Update the neighbor list with the atom positions
    nl.update(atoms)

    n_ids = nl.get_neighbors(idx)[0]

    return n_ids


def get_nn_data_from_db(dbname,n_neighbors,atom_idx,impose_neighbors=None):
    data = []
    with connect(f'hea_dbs/{dbname}') as db:

        if n_neighbors<6:
            defect_keys=['add']
        else:
            defect_keys=['none','replace']

        discards = []
        for defect_key in defect_keys:
            for row in db.select(defect=defect_key):
                atoms = db.get_atoms(row.id)

                neighbor_ids = get_nearest_neighbors(atoms,atom_idx)

                # if impose_neighbors is not None:

                #     correct_neighbors=np.all([n_id in impose_neighbors for n_id in neighbor_ids])
    
                
                
                if len(neighbor_ids)==n_neighbors and (impose_neighbors is None or np.all([n_id in impose_neighbors for n_id in neighbor_ids])):



                    if n_neighbors<6:
                        dE = db.get(slab_idx=row.slab_idx,defect='none').energy + E_bulk[row.metal] - row.energy
                    else:
                        dE = db.get(slab_idx=row.slab_idx,metal='none',defect='removed').energy + E_bulk[row.metal] - row.energy
      

                    U = dE/U_standard[row.metal]['n'] + U_standard[row.metal]['U']

                    n_symbols = atoms.get_chemical_symbols()
                    diss_metal = np.zeros(5)
                    diss_metal[metals.index(row.metal)] = 1

                    n_metals = np.zeros(5)
                    for i in neighbor_ids:
                        n_metals[metals.index(n_symbols[i])] += 1

                    datalist = np.concatenate((diss_metal,n_metals,[n_neighbors,dE,U]))

                    data.append(datalist)
                else:
                    discards.append(row.id)
                    # print(row.id,neighbor_ids)
    
    if len(discards)==0:
        print(f'No dicards in {dbname}')
    else:
        print(f'Discarded the following rows in {dbname}:', discards)
    return np.array(data)



# data fra T111
cn3 = get_nn_data_from_db('T111_HEA_out.db',3,45,impose_neighbors=[40,41,43])
cn9 = get_nn_data_from_db('T111_HEA_out.db',9,40)

# T100
cn4 = get_nn_data_from_db('T100_HEA_out.db',4,45)
cn8 = get_nn_data_from_db('T100_HEA_out.db',8,40)

# Edge
cn5 = get_nn_data_from_db('Edge_HEA_out.db',5,45)
cn7 = get_nn_data_from_db('Edge_HEA_out.db',7,1)

# Kink
cn6 = get_nn_data_from_db('Kink_HEA_out.db',6,2)


data = np.vstack((cn3,cn4,cn5,cn6,cn7,cn8,cn9))

fmt = ['%i']*11 + ['%1.6f']*2
np.savetxt('hea_data.csv',data,delimiter=',',header='Ir,Pd,Pt,Rh,Ru,n_Ir,n_Pd,n_Pt,n_Rh,n_Ru,cn,dE,U',fmt=fmt)

cns = [3,4,5,6,7,8,9]

import matplotlib
matplotlib.rcParams.update({'font.size': 16})
fig, axes = plt.subplots(ncols=5,nrows=2,figsize=(15,6),sharex=True)
# ((ax11,ax12,ax13,ax14,ax15),(ax21,ax22,ax23,ax24,ax25))
for metal, (ax1,ax2) in zip(metals,list(zip(*axes))):


    metal_mask = data[:,metals.index(metal)] == 1
    data_metal = data[metal_mask]

    for cn in cns:
        cn_mask = data_metal[:,-3] == cn
        ax_violin(ax1,ax2,data_metal[cn_mask,-2],data_metal[cn_mask,-1],cn,metal)
 


for ax in axes[0]:
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

for ax in axes[1]:
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))

    # ax.set_xlabel('Coordination Number')
    ax.set_xlabel('CN')

axes[1,3].yaxis.set_major_locator(MultipleLocator(1))


axes[0,1].set_xticks([3,4,5,6,7,8,9])
axes[0,0].set_ylabel('$\Delta E$ [eV]')
axes[1,0].set_ylabel('U [V]')


handles = []
labels = []
for metal in metals:
    handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12)) 
    labels.append('')
    handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12)) 
    labels.append('')
    handles.append(Line2D([0], [0], marker='o', color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=12)) 
    labels.append(metal)

handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12))
labels.append('')

plt.tight_layout()
plt.subplots_adjust(top=0.90)

pos1 = axes[0,0].get_position()
pos2 = axes[0,4].get_position()

# mean_pos1 = (pos1.x0 + pos1.x1)/2
# mean_pos2 = (pos2.x0 + pos2.x1)/2

fig.legend(handles=handles[1:], labels=labels[1:],
           loc='outside upper center', ncol=len(handles), mode='expand',fontsize=16,bbox_to_anchor=(pos1.x0, .5, pos2.x1-pos1.x0, 0.5),fancybox=False)


plt.savefig('DFT_cleaned_violin.png',dpi=600,bbox_inches='tight')


