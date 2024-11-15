import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import numpy as np
from utilities.compositionspace_functions import make_ternary_plot, get_molar_fractions, molar_fractions_to_cartesians


data = np.loadtxt('grid_search/grid_data.csv',delimiter=',',usecols=(0,1,2,3,4,5,6,7,32))

mfs = data[:,:-1]
Sd = data[:,-1]

# make psudo ternary
pseudo_mfs = np.hstack((np.sum(mfs[:,[1,3,5]],axis=1).reshape(-1,1),np.sum(mfs[:,[0,4]],axis=1).reshape(-1,1),np.sum(mfs[:,[2,6,7]],axis=1).reshape(-1,1)))

ternary_mfs = get_molar_fractions(0.0625,3)

pseudo_Sd = []
for mf in ternary_mfs:
    mask = np.all(np.isclose(pseudo_mfs,mf),axis=1)
    
    pseudo_Sd.append(np.max(Sd[mask]))


grid = molar_fractions_to_cartesians(ternary_mfs).T
make_ternary_plot(grid,pseudo_Sd,elements=['AuIrPt','AgPd','CuRhRu'],colormap='coolwarm_r',minval=0.0,contourlines=True)

plt.savefig('stability_plots/stability_space.png',bbox_inches='tight',dpi=600)




# Au-Cu-Pt
mask = np.isclose(np.sum(mfs[:,[1,2,5]],axis=1),1)

y=Sd[mask]
make_ternary_plot(grid,y,elements=['Au','Cu','Pt'],colormap='coolwarm_r',minval=0.0,contourlines=True,contour_levels=15)

plt.savefig('stability_plots/AuCuPt_Sd.png',bbox_inches='tight',dpi=600)

# Au-Pt edge
AuPt_mask = ternary_mfs[:,1]==0

fig,ax = plt.subplots(figsize=(3,2))

ax.plot(ternary_mfs[AuPt_mask,2],y[AuPt_mask],c='k',alpha=0.8)
ax.scatter(ternary_mfs[AuPt_mask,2],y[AuPt_mask],c='k')

ax.set_ylabel('$S_d$')
ax.set_xlabel('Au$_{1-x}$Pt$_x$')
plt.savefig('stability_plots/AuPt_Sd.png',dpi=600,bbox_inches='tight')




# Pd-Pt-Ru
mask = np.isclose(np.sum(mfs[:,[4,5,7]],axis=1),1)
y=Sd[mask]
make_ternary_plot(grid,y,elements=['Pd','Pt','Ru'],colormap='coolwarm_r',minval=0.0,contourlines=True,contour_levels=15)

plt.savefig('stability_plots/PdPtRu_Sd.png',bbox_inches='tight',dpi=600)



# Coloarbar
fig,ax = plt.subplots(figsize=(0.2,3))
cbar=plt.colorbar(ScalarMappable(cmap='coolwarm_r'),cax=ax,orientation='vertical')
cbar.set_label(label='$S_d$',rotation='horizontal')
plt.savefig('stability_plots/Sd_colorbar.png',dpi=600,bbox_inches='tight')