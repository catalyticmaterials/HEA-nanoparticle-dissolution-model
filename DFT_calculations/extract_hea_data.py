from ase.db import connect
from ase.neighborlist import NeighborList
import numpy as np
from utilities import metals


with connect('metal_dbs/bulks.db') as db:
    E_bulk={row.metal:row.energy for row in db.select()}



tol = 0.5


# Number of different slab compositions
N = 100
n_metals = len(metals)

fmt = ['%i']*n_metals*2 + ['%1.6f','%i','%i']
dataheader = ','.join(metals) + ',n_' + ',n_'.join(metals) + ',dE,CN,db_rowidx'





def get_nearest_neighbors(atoms,idx,a):

    cutoff = [a/2]*len(atoms)
    nl = NeighborList(cutoff,self_interaction=False, bothways=True)

    # Update the neighbor list with the atom positions
    nl.update(atoms)

    n_ids = nl.get_neighbors(idx)[0]

    return n_ids





def check_surface(row,surface):

    # Check if the surface is deviating from the modeled surface by restricting atoms (besides adatoms) 
    # to not have moved further than tol*a, where a is initial interatomic distance.
    atoms = row.toatoms()

    db_orig=connect(f'hea_dbs/{surface}_HEA.db')

    atoms_orig = db_orig.get_atoms(row.id)

    dists = np.linalg.norm(atoms.get_positions() - atoms_orig.get_positions(),axis=1)

    if surface=='Edge':
        a = atoms.get_cell()[1,1]/3
        
    elif surface=='Kink':
        a = atoms.get_cell()[1,1]/2.5
        
    else:
        a = atoms.get_cell()[0,0]/3


    if np.any(dists>(tol*a)) and row.defect!='adatom':
        return False



    elif np.any(dists[:-1]>(tol*a)) and row.defect=='adatom':
        return False
    

    return True


T111_fcc = np.sort([[36,37,39],[37,38,40],[38,36,41],[39,40,42],[40,41,43],[41,39,44],[42,43,36],[43,44,37],[44,42,38]],axis=1)
edge_adatom = np.array([[0,1,4,6,7],[1,2,5,7,8],[0,2,3,6,8]])
T100_hollow = np.sort([[36,37,39,40],[37,38,40,41],[38,36,41,39],[39,40,42,43],[40,41,43,44],[41,39,44,42],[42,43,36,37],[43,44,37,38],[44,42,38,36]],axis=1)

adatom_neighbors = {'T111': T111_fcc,'T100':T100_hollow,'Edge':edge_adatom}


def get_nn_data_from_db(surface,n_neighbors,target_idx):
    data = []
    discards = []
    with connect(f'hea_dbs/{surface}_HEA_out.db') as db:

        if n_neighbors<6:
            defect_keys=['adatom']
        else:
            defect_keys=['none','replace']

        
        for defect_key in defect_keys:
            for row in db.select(defect=defect_key):
                atoms = db.get_atoms(row.id)



                if surface=='Edge':
                    a = atoms.get_cell()[1,1]/3
                    
                elif surface=='Kink':
                    a = atoms.get_cell()[1,1]/2.5
                    
                else:
                    a = atoms.get_cell()[0,0]/3


                
                # Check that the target atom has the correct CN
                neighbor_ids = get_nearest_neighbors(atoms,target_idx,a)
                CN = len(neighbor_ids)
                CN_bool = CN==n_neighbors

                if n_neighbors==3:
                    pos = atoms.get_positions()
                    avec = atoms.get_cell()[0]
                    bvec = atoms.get_cell()[1]
                    min_dists=np.min([np.linalg.norm(np.delete(pos,target_idx,axis=0) - pos[target_idx] + vec,axis=1) for vec in [np.zeros(3),avec,-avec,bvec,-bvec]],axis=0)
                    neighbor_ids = np.argsort(min_dists)[:n_neighbors]
                    
                    CN_bool = np.any(np.all(np.sort(neighbor_ids).reshape(1,-1)==T111_fcc,axis=1))


                if n_neighbors<6:
                    row_diss = db.get(slab_idx=row.slab_idx,defect='none')
                else:
                    row_diss = db.get(slab_idx=row.slab_idx,metal='none',defect='vacancy')
                    
                dE = row_diss.energy + E_bulk[row.metal] - row.energy

                if CN_bool:
                    keep_bool = check_surface(row,surface)*check_surface(row_diss,surface)
                else:
                    keep_bool = False


                symbols = atoms.get_chemical_symbols()
                diss_metal = [0]*n_metals
                diss_metal[metals.index(row.metal)] = 1

                
                neighbor_metals = [0]*n_metals
                for i in neighbor_ids:
                    neighbor_metals[metals.index(symbols[i])] += 1

                datalist = diss_metal + neighbor_metals + [dE,CN,row.id-1]


                if keep_bool:
                    data.append(datalist)
                else:
                    discards.append(datalist)
    
    data = np.array(data)
    discards = np.array(discards)

    # sort data by datarow
    data = data[np.argsort(data[:,-1])]
    discards = discards[np.argsort(discards[:,-1])]
    tot = len(data)+len(discards)
    print(f'{surface},CN={n_neighbors}: kept {len(data)} ({len(data)/tot*100:1.0f}%), discarded {len(discards)} ({len(discards)/tot*100:1.0f}%)')
    np.savetxt(f'hea_dbs/discards/{surface}_CN{n_neighbors}_disc.csv',discards,delimiter=',',header=dataheader,fmt=fmt)
    return data



# data from T111
cn3 = get_nn_data_from_db('T111',3,45)
cn9 = get_nn_data_from_db('T111',9,40)

# T100
cn4 = get_nn_data_from_db('T100',4,45)
cn8 = get_nn_data_from_db('T100',8,40)

# Edge
cn5 = get_nn_data_from_db('Edge',5,45)
cn7 = get_nn_data_from_db('Edge',7,1)

# Kink
cn6 = get_nn_data_from_db('Kink',6,2)


data = np.vstack((cn3,cn4,cn5,cn6,cn7,cn8,cn9))
np.savetxt('hea_data.csv',data,delimiter=',',header=dataheader,fmt=fmt)


