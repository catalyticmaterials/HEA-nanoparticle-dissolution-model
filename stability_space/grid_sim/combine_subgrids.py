import numpy as np



subgrids = [np.loadtxt(f'grid_sim_{str(i).zfill(2)}.csv',delimiter=',') for i in range(100)]

grids_comb = np.vstack(subgrids)

fmt = ['%1.4f']*8 + ['%1.6f']*17
np.savetxt('grids_full.csv',grids_comb,delimiter=',',fmt=fmt)


