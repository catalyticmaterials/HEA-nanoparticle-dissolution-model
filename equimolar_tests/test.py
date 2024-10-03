from utilities import all_metals as metals
from utilities.particle_dissolver import Dissolver
import numpy as np 
D = Dissolver(metals,c_metals=1e-6)

for metal in metals:
    print(metal)
    print(np.around(D.U_dict[metal]['U'],decimals=4))
    print()