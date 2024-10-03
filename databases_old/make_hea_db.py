from ase.db import connect
from ase.build import fcc111,fcc100, fcc211, surface, add_adsorbate
from ase.visualize import view
import numpy as np
from copy import deepcopy
from ase.constraints import FixAtoms

lattice_parameters = np.loadtxt('lattice_parameters.csv',delimiter=',',skiprows=1)
metals = ['Ir','Pd','Pt','Rh','Ru']
lattice_parameters = {metal:a for metal,a in zip(metals,lattice_parameters)}

# Number of compositions
N=50

np.random.seed(42)
fs = np.random.dirichlet(alpha=np.ones(5), size=50)

size=(3,3,5)


for idx,f in enumerate(fs):
    # print(f)
    symbols = np.random.choice(metals, size=np.product(size), p=f, replace=True)
    a = np.mean([lattice_parameters[symbol] for symbol in symbols])
    
    
    with connect('hea_dbs/T111_HEA.db') as db:
        T111 = fcc111('X',size,a,vacuum=10)
        T111.set_chemical_symbols(symbols)
        # Fix all but the two top layers of atoms
        constraint = FixAtoms(indices=range(27))
        T111.set_constraint(constraint)
        
        # Metal of interest at atom index 40
        moi = symbols[40]

        db.write(T111,slab_idx=idx,metal=moi,defect='none')
    
        symbols_copy = deepcopy(symbols)



        # Remove
        T111_remove = deepcopy(T111)
        del T111_remove[40]
        db.write(T111_remove,slab_idx=idx,metal='none',defect='removed')

        for metal in metals:
            
            # Replace center metal with metal
            if moi!=metal:
                T111_replace = deepcopy(T111)
                symbols_copy[40] = metal
                T111_replace.set_chemical_symbols(symbols_copy)
                db.write(T111_replace,slab_idx=idx,metal=metal,defect='replace')

                

            # Add
            T111_add = deepcopy(T111)
            add_adsorbate(T111_add,metal,position='fcc',offset=(1,1),height=a/np.sqrt(3))
            db.write(T111_add,slab_idx=idx,metal=metal,defect='add')

            
            


    

    with connect('hea_dbs/T100_HEA.db') as db:

        T100 = fcc100('X',(3,3,5),a=a,vacuum=10)
        T100.set_chemical_symbols(symbols)
        constraint = FixAtoms(indices=range(27))
        T100.set_constraint(constraint)
        
   
        
        # Metal of interest at atom index 40
        moi = symbols[40]

        db.write(T100,slab_idx=idx,metal=moi,defect='none')
    
        symbols_copy = deepcopy(symbols)

        # Remove
        T100_remove = deepcopy(T100)
        del T100_remove[40]
        db.write(T100_remove,slab_idx=idx,metal='none',defect='removed')

        
        for metal in metals:
            
            # Replace center metal with metal
            if moi!=metal:
                T100_replace = deepcopy(T100)
                symbols_copy[40] = metal
                T100_replace.set_chemical_symbols(symbols_copy)
                db.write(T100_replace,slab_idx=idx,metal=metal,defect='replace')

                

            # Add
            T100_add = deepcopy(T100)
            add_adsorbate(T100_add,metal,position='hollow',offset=(1,1),height=a/2)
            db.write(T100_add,slab_idx=idx,metal=metal,defect='add')



    with connect('hea_dbs/Edge_HEA.db') as db:

        # Edge
        E = fcc211('X',size,a=a,vacuum=10)
        E.set_chemical_symbols(symbols)

        constraint = FixAtoms(indices=range(18,45))
        E.set_constraint(constraint)

        # Metal of interest at atom index 1
        moi = symbols[1]

        symbols_copy = deepcopy(symbols)

        db.write(E,slab_idx=idx,metal=moi,defect='none')

        # Remove
        E_remove = deepcopy(E)
        del E_remove[1]
        db.write(E_remove,slab_idx=idx,metal='none',defect='removed')

        for metal in metals:
            
            # Replace center metal with metal
            if moi!=metal:
                E_replace = deepcopy(E)
                symbols_copy[1] = metal
                E_replace.set_chemical_symbols(symbols_copy)
                db.write(E_replace,slab_idx=idx,metal=metal,defect='replace')
                

                

            # Add
            E_add = deepcopy(E)
            add_adsorbate(E_add,metal,height=a*np.sqrt(6)/12,position=(0,0),offset=(0,1/3))
            db.write(E_add,slab_idx=idx,metal=metal,defect='add')             


    with connect('hea_dbs/Kink_HEA.db') as db:

        # Kink
        K = fcc211('X',size,a=a,vacuum=10)
        K.set_chemical_symbols(symbols)
        
        cell = K.get_cell()

        a_,b,c = cell
        
        b[0] = -a_[0]/3
        b[1] = K.positions[2][1]
        h = K.positions[3] - K.positions[6]
        b[2] = -h[2]

        # view(K)
        K.set_cell(np.array([a_,b,c]))
        constraint = FixAtoms(indices=range(18,45))
        K.set_constraint(constraint)

        del K[[8,14,26,32,20,38,44]]

        symbols = K.get_chemical_symbols()
        
        # Metal of interest at atom index 2
        moi = symbols[2]

        symbols_copy = deepcopy(symbols)

        db.write(K,slab_idx=idx,metal=moi,defect='none')

        # Remove
        K_remove = deepcopy(K)
        del K_remove[2]
        db.write(K_remove,slab_idx=idx,metal='none',defect='removed')


        for metal in metals:
            
            # Replace center metal with metal
            if moi!=metal:
                K_replace = deepcopy(K)
                symbols_copy[2] = metal
                K_replace.set_chemical_symbols(symbols_copy)
                db.write(K_replace,slab_idx=idx,metal=metal,defect='replace')