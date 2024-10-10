import numpy as np
from utilities import metals


data  = np.loadtxt('grid_sim/grids_full.csv',delimiter=',')
data_isc = np.loadtxt('grid_data_.csv',delimiter=',',usecols=(8,9,10,11,12,13,14,15))

newdata = np.hstack((data[:,:8],data_isc,data[:,8:]))

header = ','.join(metals) + ',' + '_isc,'.join(metals) + '_isc,' + '_fsc,'.join(metals) + '_fsc,' + '_d,'.join(metals) + '_d,Sd'
np.savetxt('grid_data.csv',newdata,delimiter=',',fmt='%1.6f',header=header)