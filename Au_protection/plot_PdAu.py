import numpy as np
import matplotlib.pyplot as plt
from utilities.colors import metal_colors

Au,Pd,Au_f,Pd_f,Sd,Sd_se = np.loadtxt('PdAu.csv',delimiter=',',skiprows=1).T


fig,ax = plt.subplots(figsize=(4,3))

ax2 = ax.twinx()
ax2.errorbar(Au,Sd,Sd_se,capsize=1.2,c='k',capthick=0.5,elinewidth=0.8,linewidth=1)
ax.bar(Au,Pd_f*100,color=metal_colors['Pd'],width=0.045,label='Pd')
ax.bar(Au,Au_f*100,bottom=Pd_f*100,color=metal_colors['Au'],width=0.045,label='Au')

ax.set_xlabel('$Pd_xAu_{1-x}$')
ax.set_ylabel('Composition (at. %)')
ax2.set_ylabel('$S_d$')
ax.set_xlim(-0.025,1.025)
ax.set_ylim(0,100)
ax2.set_ylim(0,1)
plt.tight_layout()

plt.savefig('PdAu.png',dpi=600,bbox_inches='tight')