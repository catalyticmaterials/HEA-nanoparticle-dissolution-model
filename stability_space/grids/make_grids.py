from utilities.compositionspace_functions import get_molar_fractions
import numpy as np


grid = get_molar_fractions(0.0625,8)

grid_split = np.array_split(grid,100)


for i,sub_grid in enumerate(grid_split):
    idx = str(i).zfill(2)
    np.savetxt(f'grids/subgrid_{idx}.csv',sub_grid,delimiter=',',fmt='%1.4f')