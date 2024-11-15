from utilities.colors import metal_colors
from utilities import metals
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator

n_metals = len(metals)

data = np.loadtxt('shape/shape_dependence.csv',delimiter=',',skiprows=1)

N = data[:,0]
E = data[:,1]
comp111 = data[:,2:n_metals+2]
comp111_std = data[:,n_metals+2:2*n_metals+2]
comp_diss = data[:,2*n_metals+2:3*n_metals+2]
comp_diss_std = data[:,3*n_metals+2:-2]
df = data[:,-2]
df_std = data[:,-1]



fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)
axes = axes.flatten()
for i,metal in enumerate(metals):

    axes[i].errorbar(range(len(N)),comp111[:,i],comp111_std[:,i],fmt='.',color=metal_colors[metal],capsize=5)

    

axes[0].set_xticks(range(len(N)),labels=[1.0,1.1,1.2,1.3,1.4,1.5])
axes[0].yaxis.set_minor_locator(MultipleLocator(0.1))

for i in range(4,8):
    axes[i].set_xlabel('$E_{(100)}/E_{(111)}$')
axes[0].set_ylabel('Surface comp.')
axes[4].set_ylabel('Surface comp.')


handles = []
for metal in metals:
    handles.append(Line2D([0], [0], marker='o', color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=10)) 


plt.tight_layout()
pos1 = axes[0].get_position()
pos2 = axes[3].get_position()


fig.legend(handles=handles, labels=metals,
           loc='outside upper center', ncol=n_metals, mode='expand',bbox_to_anchor=(pos1.x0, .5, pos2.x1-pos1.x0, 0.5),fancybox=False)


# plt.minorticks_off()


plt.subplots_adjust(top=0.9)
plt.savefig('shape/shape_dependence_111_comp.png',dpi=600)





fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)
axes = axes.flatten()
for i,metal in enumerate(metals):

    axes[i].errorbar(range(len(N)),comp_diss[:,i],comp_diss_std[:,i],fmt='.',color=metal_colors[metal],capsize=5)


# axes[1,0].set_xscale('log')

axes[0].set_xticks(range(len(N)),labels=[1.0,1.1,1.2,1.3,1.4,1.5])
axes[0].yaxis.set_minor_locator(MultipleLocator(0.05))

for i in range(4,8):
    axes[i].set_xlabel('$E_{(100)}/E_{(111)}$')

axes[0].set_ylabel('Dissolved comp.')
axes[4].set_ylabel('Dissolved comp.')


handles = []
for metal in metals:
    handles.append(Line2D([0], [0], marker='o', color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=10)) 

plt.tight_layout()
pos1 = axes[0].get_position()
pos2 = axes[3].get_position()

fig.legend(handles=handles, labels=metals,
           loc='outside upper center', ncol=n_metals, mode='expand',bbox_to_anchor=(pos1.x0, .5, pos2.x1-pos1.x0, 0.5),fancybox=False)


# plt.minorticks_off()
plt.subplots_adjust(top=0.9)
plt.savefig('shape/shape_dependence_diss_comp.png',dpi=600)




rcParams.update({'font.size': 12})

fig,ax = plt.subplots(figsize=(3,3))
ax.errorbar(E,df,df_std,fmt='.')

ax.set_xlabel('$E_{(100)}/E_{(111)}$')
ax.set_ylabel('$S_d$')


ax.xaxis.set_major_locator(MultipleLocator(0.1))
# ax.yaxis.set_major_locator(MultipleLocator(0.05))



# plt.tight_layout()
plt.savefig('shape/Sd_shape_dependence.png',dpi=600,bbox_inches='tight')

