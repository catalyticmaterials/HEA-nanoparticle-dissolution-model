from ase.db import connect
from ase.build import fcc111,fcc100, fcc211, surface, add_adsorbate
from ase.constraints import FixAtoms
from ase.visualize import view
import numpy as np
from copy import deepcopy

lattice_parameters = np.loadtxt('lattice_parameters.csv',delimiter=',',skiprows=1)

with connect('metal_dbs/metals_dissolution.db') as db:
	for metal,a in zip(('Ag','Au','Cu','Ir','Pd','Pt','Rh','Ru'),lattice_parameters):

		# Terrace 111
		T1 = fcc111(metal, (3,3,5),a=a, vacuum=10)
		constraint = FixAtoms(indices=range(27))
		T1.set_constraint(constraint)
		
		T1_rm = deepcopy(T1)
		
		del T1_rm[40]

		T1_add = deepcopy(T1)

		add_adsorbate(T1_add,metal,position='fcc',offset=(1,1),height=a/np.sqrt(3))

		db.write(T1,metal=metal,surface='T111',defect='none')
		db.write(T1_rm,metal=metal,surface='T111',defect='vacancy')
		db.write(T1_add,metal=metal,surface='T111',defect='adatom')


		# Terrace 100
		T2 = fcc100(metal,(3,3,5),a=a,vacuum=10)
		T2.set_constraint(constraint)

		T2_rm = deepcopy(T2)
		
		del T2_rm[40]
		
		# Add on 100
		T2_add = deepcopy(T2)

		add_adsorbate(T2_add,metal,position='hollow',offset=(1,1),height=a/2)

		db.write(T2,metal=metal,surface='T100',defect='none')
		db.write(T2_rm,metal=metal,surface='T100',defect='vacancy')
		db.write(T2_add,metal=metal,surface='T100',defect='adatom')


		
		# Edge
		E = fcc211(metal,(3,3,5),a=a,vacuum=10)
		constraint = FixAtoms(indices=range(18,45))
		E.set_constraint(constraint)

		E_rm = deepcopy(E)
		del E_rm[1]
		
		E_add = deepcopy(E)

		add_adsorbate(E_add,metal,height=a*np.sqrt(6)/12,position=(0,0),offset=(0,1/3))

		db.write(E,metal=metal,surface='edge',defect='none')
		db.write(E_rm,metal=metal,surface='edge',defect='vacancy')
		db.write(E_add,metal=metal,surface='edge',defect='adatom')



		# Kink
		K = fcc211(metal,(3,3,5),a=a,vacuum=10)
		K.set_constraint(constraint)

		cell = K.get_cell()

		a_,b,c = cell

		b[0] = -a_[0]/3
		b[1] = K.positions[2][1]
		h = K.positions[3] - K.positions[6]
		b[2] = -h[2]

		# view(K)
		K.set_cell(np.array([a_,b,c]))
		

		del K[[8,14,26,32,20,38,44]]

		K_rm = deepcopy(K)

		del K_rm[2]

		db.write(K,metal=metal,surface='kink',defect='none')
		db.write(K_rm,metal=metal,surface='kink',defect='vacancy')

		

		