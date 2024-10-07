import numpy as np



subgrids = [np.loadtxt(f'grid_sim_{str(i).zfill(2)}.csv',delimiter=',') for i in range(100)]

grids_comb = np.vstack(grids)

np.savetxt('grids_full.csv',grids_comb,delimiter=',',fmt='%1.6f')


