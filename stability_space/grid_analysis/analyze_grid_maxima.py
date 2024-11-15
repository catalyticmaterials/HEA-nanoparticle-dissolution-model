import numpy as np
from utilities.compositionspace_functions import get_molar_fractions_around

# Sort grid maxima into regions


full_data = np.loadtxt('grid_search/grid_data.csv',delimiter=',',usecols=(0,1,2,3,4,5,6,7,32))

max_data = np.loadtxt('grid_analysis/grid_maxima.csv',delimiter=',')


mfs_max = max_data[:,:-1]
Sd_max = max_data[:,-1]

regions = []
regions_Sd = []
while len(mfs_max)>0:

    if len(mfs_max)==1:
        regions.append(np.atleast_2d(mfs_max))
        regions_Sd.append(Sd_max[0])
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


region_ids = np.concatenate([[i]*len(region) for i,region in enumerate(regions)]).reshape(-1,1)
regions_Sd_data = np.array([regions_Sd[i] for i in region_ids.ravel()]).reshape(-1,1)
regions_data = np.vstack(regions)

data = np.hstack((regions_data*100,regions_Sd_data,region_ids))
np.savetxt('grid_analysis/max_regions.csv',data,delimiter=',',fmt=['%1.1f']*8+['%1.2f','%i'],header='Ag,Au,Cu,Ir,Pd,Pt,Rh,Ru,Sd,idx')




# Inspect last region
mfs = full_data[:,:-1]
Sd = full_data[:,-1]

mfs_around = get_molar_fractions_around(regions[-1][0],0.0625)

ids = np.array([np.where(np.all(np.isclose(mfs,mf_around),axis=1))[0][0] for mf_around in mfs_around])



Sd_around = Sd[ids]

print()
sort_ids = np.argsort(-Sd_around)

print(mfs_around[sort_ids[:3]])
print(Sd_around[sort_ids[:3]])

