import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

data = np.loadtxt('grid_sim/grids_full.csv',delimiter=',',usecols=(0,2,4,6,7,24))

mfs = data[:,:-1]
Sd = data[:,-1]

mask = np.isclose(np.sum(mfs,axis=1),1)

mfs = mfs[mask]
Sd = Sd[mask]
data = data[mask]

step_size = 0.0625*2
eps = 1e-6
# step_size+=eps

n=len(mfs)


def maxima(i):
    
    if Sd[i]==0.0:
        return False
    
    delta_mf =np.sum(np.abs(mfs-mfs[i]),axis=1)/2 
    mask = np.isclose(delta_mf,step_size)
    # mask[i]=False

    Sd_around = Sd[mask]

    if np.all(Sd[i]>=Sd_around):
        sort_ids = np.argsort(-Sd_around)
        print(mfs[mask][sort_ids[:3]])
        print(Sd_around[sort_ids[:3]])

    return np.all(Sd[i]>=Sd_around)
        

if __name__ == '__main__':
    with Pool(4) as pool:
        maxima_mask = list(tqdm(pool.imap(maxima,range(n)),total=n, desc= 'Processing',mininterval=10))


    maxima_mask = np.array(maxima_mask)

    np.savetxt('grid_maxima_nn.csv',data[maxima_mask],header='Ag,Cu,Pd,Rh,Ru,Sd',fmt='%1.6f',delimiter=',')

    

    