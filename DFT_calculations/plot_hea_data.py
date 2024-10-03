import matplotlib.pyplot as plt
import numpy as np
from utilities import metals,U_standard
from utilities.colors import metal_colors
from matplotlib.ticker import MultipleLocator

n_metals = len(metals)

fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)

data = np.loadtxt('hea_data.csv',skiprows=1,delimiter=',')
target = data[:,:n_metals]

CNs = np.arange(3,10)

for i in range(4):
    axes[1,i].set_xlabel('CN')
axes[0,0].set_ylabel(r'$\Delta E [eV]$')
axes[1,0].set_ylabel(r'$\Delta E [eV]$')
axes[0,0].set_xticks(CNs,CNs)
axes[0,0].yaxis.set_minor_locator(MultipleLocator(1))

for metal,ax in zip(metals,axes.flatten()):

    metal_idx = metals.index(metal)
    
    metal_mask = target[:,metal_idx]==1
    metal_data = data[metal_mask]
    dE = metal_data[:,-3]
    CN = metal_data[:,-2]

    # print(metal)
    for cn in CNs:
        CN_mask = cn==CN
        parts = ax.violinplot(dE[CN_mask],[cn],showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(metal_colors[metal])
            pc.set_edgecolor(metal_colors[metal])
            pc.set_alpha(1)
        
        mu=np.mean(dE[CN_mask])
        s = np.std(dE[CN_mask])

        z = np.abs((dE[CN_mask]-mu)/s)
        # print(cn)
        CN_data = metal_data[CN_mask]
        # if np.any(z>3):
        #     print(metal,cn)
        #     print(CN_data[z>3])
        if cn==6:
            print(metal,np.sum(CN_data[:,8+metal_idx]>0))



plt.tight_layout()
plt.savefig('hea_data.png',dpi=600,bbox_inches='tight')