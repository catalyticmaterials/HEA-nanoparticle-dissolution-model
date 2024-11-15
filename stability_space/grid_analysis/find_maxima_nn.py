import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

# find maxima in grid exluding Au, Ir, and Pt


data = np.loadtxt('grid_search/grid_data.csv',delimiter=',',usecols=(0,2,4,6,7,32))

mfs = data[:,:-1]
Sd = data[:,-1]

mask = np.isclose(np.sum(mfs,axis=1),1)

mfs = mfs[mask]
Sd = Sd[mask]
data = data[mask]

step_size = 0.0625


n=len(mfs)


def maxima(i):
    # Check if grid point is has Sd above or equal to surrounding grid points.
    delta_mf =np.sum(np.abs(mfs-mfs[i]),axis=1)/2 
    mask = (delta_mf<=step_size)
    mask[i]=False

    Sd_around = Sd[mask]

    return np.all(Sd[i]>=Sd_around)
        

if __name__ == '__main__':
    with Pool(4) as pool:
        maxima_mask = list(tqdm(pool.imap(maxima,range(n)),total=n, desc= 'Processing',mininterval=10))


    maxima_mask = np.array(maxima_mask)

    np.savetxt('grid_analysis/grid_maxima_nn.csv',data[maxima_mask],header='Ag,Cu,Pd,Rh,Ru,Sd',fmt='%1.6f',delimiter=',')

    

    