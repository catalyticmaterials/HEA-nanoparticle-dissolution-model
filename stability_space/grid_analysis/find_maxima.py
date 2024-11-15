import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

data = np.loadtxt('../grid_search/grids_data.csv',delimiter=',',usecols=(0,1,2,3,4,5,6,7,32))

mfs = data[:,:-1]
Sd = data[:,-1]

step_size = 0.0625
eps = 1e-6
step_size+=eps

n=len(data)


def maxima(i):
    # Check if grid point is has Sd above or equal to surrounding grid points.
    delta_mf =np.sum(np.abs(mfs-mfs[i]),axis=1)/2 
    mask = (delta_mf<=step_size)
    mask[i]=False

    Sd_around = Sd[mask]

    return np.all(Sd[i]>=Sd_around)
        

if __name__ == '__main__':
    with Pool(32) as pool:
        maxima_mask = list(tqdm(pool.imap(maxima,range(n)),total=n, desc= 'Processing',mininterval=10))


    maxima_mask = np.array(maxima_mask)

    np.savetxt('grid_analyis/grid_maxima.csv',data[maxima_mask],header='Ag,Au,Cu,Ir,Pd,Pt,Rh,Ru,Sd',fmt='%1.6f',delimiter=',')

    
