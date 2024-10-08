import numpy as np

full_data = np.loadtxt('grid_sim/grids_full.csv',delimiter=',',usecols=(0,1,2,3,4,5,6,7,24))

max_data = np.loadtxt('grid_maxima.csv',delimiter=',')


mfs_max = max_data[:,:-1]
Sd_max = max_data[:,-1]

regions = []
regions_Sd = []
while len(mfs_max)>0:

    if len(mfs_max)==1:
        regions.append(np.atleast_2d(mfs_max))
        regions_Sd.append(Sd_max)
        break

    if np.any(mfs_max==1):
        
        idx = np.where(np.any(mfs_max==1,axis=1))[0][0]
    
    else:
        idx=0

    current_region = np.array([mfs_max[idx]])
    mfs_max = np.delete(mfs_max,idx,axis=0)
    regions_Sd.append(Sd_max[idx])
    Sd_max = np.delete(Sd_max,idx)

    while True:
        add_ids = np.where([np.any(np.sum(np.abs(current_region-mf),axis=1)/2==0.0625) for mf in mfs_max])[0]
        if len(add_ids)==0:
            break

        current_region = np.vstack((current_region,np.atleast_2d(mfs_max[add_ids])))
        mfs_max = np.delete(mfs_max,add_ids,axis=0)
        Sd_max = np.delete(Sd_max,add_ids)

    regions.append(current_region)

print(regions)
print(regions_Sd)

from utilities.compositionspace_functions import get_molar_fractions_around

mfs = full_data[:,:-1]
Sd = full_data[:,-1]

mfs_around = get_molar_fractions_around(regions[-1][0],0.0625)

ids = np.array([np.where(np.all(np.isclose(mfs,mf_around),axis=1))[0][0] for mf_around in mfs_around])



Sd_around = Sd[ids]

print()
sort_ids = np.argsort(-Sd_around)

print(mfs_around[sort_ids[:3]])
print(Sd_around[sort_ids[:3]])



import matplotlib.pyplot as plt

mask = np.isclose(np.sum(mfs[:,[1,5]],axis=1),1)

plt.errorbar(mfs[mask,1],Sd[mask],0.01)

plt.show()