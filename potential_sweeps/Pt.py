import numpy as np
from utilities.particle_dissolver import Dissolver
from utilities import metals,load_regressor
from utilities.compositionspace_functions import get_molar_fractions
import matplotlib.pyplot as plt
from ase.visualize import view

Ns= [25,50,100,201,405,711,1289,1925,3355] # approximate diameters (nm): 0.7,1.1, 1.2, 1.8, 2.3, 2.8, 3.5, 4.0, 5.0:

model_path = '../regression_model/IrPdPtRhRu_multilinear.regressor'
regressor = load_regressor(model_path)

dissolver = Dissolver(c_metals=1e-6)
dissolver.set_regressor(regressor)
# dissolver.set_feature_generator('LinearNeighborID')

U_step = 0.1
U_fine_step = 0.01
potentials=np.arange(-0.2,2.0,U_step)



def dissolve(N):

    comp = np.array([0.0,0.0,1.0,0.0,0.0])
    np.random.seed(42)
    particle = dissolver.make_particle(comp,N,return_particle=True)

   
    # First take big potential steps
    for U in potentials:
        
        atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
        N_final_111 = len(dissolver.get_111_atoms(atoms))

        if N_final_111==0.0:
            # Go back to previous step and break
            U_min = U-U_step
            break

    U = U_min + U_fine_step
    # Increment U in smaller steps
    while True:
        atoms,_ = dissolver.dissolve_atoms(U,atoms=particle.copy(),relax_cn=True)
        N_final_111 = len(dissolver.get_111_atoms(atoms))

        if N_final_111==0:
            return U
        else:
            U+=U_fine_step




Uc = [dissolve(N) for N in Ns]
print(Uc)
D = np.array([0.7,1.1,1.2,1.8, 2.3, 2.8, 3.5, 4.0, 5.0])

r_ = 2/(D/2)

plt.scatter(r_,Uc)

plt.xlabel('2/r')
plt.ylabel('$U_{Pt}^0$ [V]')

plt.show()