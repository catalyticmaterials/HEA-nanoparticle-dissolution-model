from utilities.compositionspace_functions import make_ternary_plot, molar_fractions_to_cartesians
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from utilities import U_diss

data = np.loadtxt('AuCuPt_potential_sweep.csv',delimiter=',',skiprows=1)

metals = ['Au','Cu','Pt']

grid = molar_fractions_to_cartesians(data[:,:3]).T
Ucs = data[:,3:]


U=np.max(Ucs,axis=1)



make_ternary_plot(grid,U,metals,vmin=np.min(U),vmax=np.max(U),colorbar=True,contour_levels=30,cbar_label='$U_c$ [V]',colormap='magma',maxval=0.95,minval=0.0,cbar_ticks=np.arange(0.2,1.41,0.2))
print(data[np.argmax(U),:3],np.max(U))
plt.savefig('AuCuPt_Ucrit.png',dpi=600,bbox_inches='tight')

from utilities import load_regressor

reg = load_regressor('../utilities/AgAuCuIrPdPtRhRu_multilinear.regressor')

U_diss_pred = {metal:np.round(U_diss(0.0,metal,1e-6),decimals=3) for metal in metals}
print(U_diss_pred)


stop


potentials=np.arange(0.2,1.1,0.01)

fd_data = np.loadtxt('AuCuPt_potential_sweep_fd.csv',delimiter=',',skiprows=1)

grid = molar_fractions_to_cartesians(fd_data[:,:3]).T
fd_all = fd_data[:,3:]


sm = ScalarMappable(cmap='coolwarm_r')

fig,axes = plt.subplots(ncols=3,figsize=(9,3))

for metal,ax in zip(metals,axes):
    U_above = potentials[potentials >= U_diss_pred[metal]]
    print(metal,U_above[0])
    idx = np.where(potentials==U_above[0])[0][0]
    
    fd = fd_all[:,idx]

    make_ternary_plot(grid,fd,metals,ax=ax,vmin=0.0,vmax=1.0,colormap='coolwarm_r',minval=0.0)


plt.subplots_adjust(bottom=0.1)
cbar_ax = fig.add_axes([0.2, 0.0, 0.6, 0.06])

cbar = plt.colorbar(sm,cax=cbar_ax,ticks=np.arange(0,1.1,0.2),orientation='horizontal')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$S_d$',size=14)
plt.tight_layout()

plt.savefig(f'AuCuPt_Ucrit_fd.png',dpi=600,bbox_inches='tight')




U_max = np.max(U)
print(U_max)
idx = np.where(np.isclose(potentials,U_max))[0][0]
fd = fd_all[:,idx]

fig,ax=make_ternary_plot(grid,fd,metals,vmin=0.0,vmax=1.0,colormap='coolwarm_r',minval=0.0)
plt.colorbar(sm,ax=ax,shrink=0.5,anchor=(0.0,0.85),ticks=np.arange(0,1.1,0.2),label='$S_d$')

plt.savefig(f'AuCuPt_max_Ucrit_fd.png',dpi=600,bbox_inches='tight')