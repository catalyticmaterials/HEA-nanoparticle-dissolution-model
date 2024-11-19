from ase.db import connect
import numpy as np
from ase.neighborlist import NeighborList, natural_cutoffs


# Script to check for deviations in surface structures due to the structure relatation
# Ignores adatom, they will be allowed to move more than the surface, but will be checked by CN and fcc position for (111)
# A surface will be deemed to deviate from modelled surface and discarded as an outlier, if any atom has moved further than an half the initials interatomic distance. 


tol = 0.5

for surface in ('T111','T100','Edge','Kink'):
    
    
    print(surface)
    with connect(f'hea_dbs/{surface}_HEA.db') as db_orig, connect(f'hea_dbs/{surface}_HEA_out.db') as db:#,connect(f'hea_dbs/{surface}_HEA_check.db') as db_check, connect('hea_dbs/discarded_structures.db') as db_discard:
        count=0
        for row in db.select():
            

            atoms_orig = db_orig.get_atoms(row.id)
            atoms = db.get_atoms(row.id)

            dists = np.linalg.norm(atoms.get_positions() - atoms_orig.get_positions(),axis=1)

            if surface=='Edge':
                a = atoms.get_cell()[1,1]/3
                continue
                
            elif surface=='Kink':
                a = atoms.get_cell()[1,1]/2.5
                
            else:
                a = atoms.get_cell()[0,0]/3
                continue

            # row_kwargs = {key:row[key] for key in row._keys}

            if np.any(dists>(tol*a)) and row.defect!='adatom' and (surface!='Kink' or row.defect=='vacancy'):
                mask = dists>tol*a
                print(row.id,row.slab_idx,row.defect,np.where(mask),np.around(dists[mask]/a,decimals=2))
                # view(atoms)
                count+=1

            elif np.any(dists[:-1]>(tol*a)) and row.defect=='adatom':
                mask = dists[:-1]>tol*a
                print(row.id,row.slab_idx,row.defect, np.where(mask),np.around(dists[:-1][mask]/a,decimals=2))
                # view(atoms)
                count+=1
            
            # Treat kink atom as adatom
            elif surface=='Kink' and np.any(np.delete(dists,2)>tol*a) and row.defect!='vacancy':
                dists_ = np.delete(dists,2)
                mask = dists_>tol*a
                print(row.id,row.slab_idx,row.defect, np.where(mask),np.around(dists_[mask]/a,decimals=2))
                count+=1


            if row.defect!='vacancy':
                cutoff = natural_cutoffs(atoms)


                nl = NeighborList(cutoff,self_interaction=False, bothways=True)

                # Update the neighbor list with the atom positions
                nl.update(atoms)

                n_ids = nl.get_neighbors(2)[0]
                if len(n_ids)==6:
                    if np.any(np.sort(n_ids)!=np.array([0,1,3,5,6,10])):
                        print(row.id,n_ids)
            

    print('count:',count)



    print()