import sys
sys.path.append('..')

from ase.db import connect
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scripts import metals, U_standard
from scripts.colors import metal_colors



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



# def get_dE_and_U(db_path,adatom=True):
#     with connect(db_path) as db: 
        
#         dE_rem = np.array([db.get(slab_idx=i,metal='none',defect='removed').energy + E_bulk[metal] - util_func(db,i,metal) for i in range(N)])
#         if adatom:
#             dE_ad = np.array([db.get(slab_idx=i,defect='none').energy + E_bulk[metal] - db.get(slab_idx=i,metal=metal,defect='add').energy for i in range(N)])

#     U_rem = dE_rem/U_metal['n'] + U_metal['U']

#     if adatom:
#         U_ad = dE_ad/U_metal['n'] + U_metal['U']
#         return dE_rem, U_rem, dE_ad, U_ad
    
#     return dE_rem, U_rem
    



def get_nearest_neighbors(atoms,idx,mult=0.9):
    cutoff = natural_cutoffs(atoms,mult=mult)


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


                    neighbor_features = [0]*5*6
                    for i in neighbor_ids:

                        next_neighbors = len(get_nearest_neighbors(atoms,i,mult=1.0))
                        
                        neighbor_features[(next_neighbors-7)*5 + metals.index(n_symbols[i])] += 1

                    datalist = np.concatenate((diss_metal,neighbor_features,[n_neighbors,dE,U]))

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

fmt = ['%i']*36 + ['%1.6f']*2
np.savetxt('hea_data_ncn.csv',data,delimiter=',',header='Ir,Pd,Pt,Rh,Ru,neighbor_cn7,neighbor_cn8,neighbor_cn9,neighbor_cn10,neighbor_cn11,neighbor_cn12,cn,dE,U',fmt=fmt)
