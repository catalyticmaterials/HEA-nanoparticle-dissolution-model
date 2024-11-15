import numpy as np
import matplotlib.pyplot as plt
from utilities import metals, U_diss
from utilities.colors import metal_colors
from matplotlib.ticker import MultipleLocator



U,fd,fd_se,std = np.loadtxt(f'potential/potential_dependence.csv',delimiter=',',skiprows=1).T

fig,ax = plt.subplots(figsize=(6,3.5))

ax.errorbar(U,fd,std,fmt='.',label='$S_d$',markersize=4,c='orangered',capsize=3,elinewidth=0.5,zorder=1,ecolor='k',capthick=1,markeredgecolor='k',mew=0.5)

for metal in metals:
    ax.axvline(U_diss(0.0,metal,1e-6),color=metal_colors[metal],label=metal,zorder=0)

ax.set_xlabel('U [V]')
ax.set_ylabel('$S_d$')
ax.legend(loc='lower left')

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

plt.savefig(f'potential/potential_dependence.png',dpi=600,bbox_inches='tight')
