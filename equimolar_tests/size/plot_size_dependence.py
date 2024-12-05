from utilities.colors import metal_colors
from utilities import metals
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams


n_metals = len(metals)

data = np.loadtxt('size/size_dependence.csv',delimiter=',',skiprows=1)

N = data[:,0]
comp111 = data[:,1:n_metals+1]
comp111_std = data[:,n_metals+1:2*n_metals+1]
comp_diss = data[:,2*n_metals+1:3*n_metals+1]
comp_diss_std = data[:,3*n_metals+1:-2]
df = data[:,-2]
df_std = data[:,-1]


fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)
axes = axes.flatten()
for i,metal in enumerate(metals):

    axes[i].errorbar(range(len(N)),comp111[:,i],comp111_std[:,i],fmt='.',color=metal_colors[metal],capsize=5)

    

axes[0].set_xticks(range(len(N)),labels=[1.8, 2.3, 2.8, 3.5, 4.0, 5.0])
axes[0].yaxis.set_minor_locator(MultipleLocator(0.1))

for i in range(4,8):
    axes[i].set_xlabel('Particle size [nm]')
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




plt.subplots_adjust(top=0.9)
plt.savefig('size/size_dependence_111_comp.png',dpi=600)





fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)
axes = axes.flatten()
for i,metal in enumerate(metals):

    axes[i].errorbar(range(len(N)),comp_diss[:,i],comp_diss_std[:,i],fmt='.',color=metal_colors[metal],capsize=5)




axes[0].set_xticks(range(len(N)),labels=[1.8, 2.3, 2.8, 3.5, 4.0, 5.0])
axes[0].yaxis.set_minor_locator(MultipleLocator(0.05))

for i in range(4,8):
    axes[i].set_xlabel('Particle size [nm]')

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



plt.subplots_adjust(top=0.9)
plt.savefig('size/size_dependence_diss_comp.png',dpi=600)








rcParams.update({'font.size': 12})

fig,ax = plt.subplots(figsize=(3,3))
ax.errorbar(N,df,df_std,fmt='.')

ax.set_xlabel('Particle size (nm)')
ax.set_ylabel('$S_d$')
ax.set_xscale('log')
ax.set_xticks(N,labels=[1.8, 2.3, 2.8, 3.5, 4.0, 5.0])

plt.minorticks_off()
plt.savefig('size/Sd_size_dependence.png',dpi=600,bbox_inches='tight')



max_error = 0.05
max_sd_error = 0.01

fig,(ax,ax3) = plt.subplots(nrows=2,figsize=(4,5),sharex=True)

norm_error_comp111 = np.linalg.norm(comp111_std,axis=1)
norm_error_comp_diss = np.linalg.norm(comp_diss_std,axis=1)

ax.scatter(N,norm_error_comp111,label=r'$\sigma$(111)')
ax.scatter(N,norm_error_comp_diss,label=r'$\sigma$(diss.)')

ax2 = ax.twinx()
n_samples_111 = np.ceil((norm_error_comp111/max_error)**2)
n_samples_diss = np.ceil((norm_error_comp_diss/max_error)**2)
ax2.plot(N,n_samples_111,label='$N_{samples}$(111)')
ax2.plot(N,n_samples_diss,label='$N_{samples}$(diss.)')


h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2)

ax3.set_xlabel('Particle size (nm)')
ax.set_ylabel(r'$\sigma_{composition}$')
ax.set_xscale('log')
ax.set_xticks(N,labels=[1.8, 2.3, 2.8, 3.5, 4.0, 5.0])



ax3.set_ylabel(r'$\sigma_{S_d}$')
ax3.scatter(N,df_std,label=r'$\sigma_{S_d}$',color='tab:green')
ax4 = ax3.twinx()
n_samples_df = np.ceil((df_std/max_sd_error)**2)
ax4.plot(N,n_samples_df,label='$N_{samples}$',color='tab:green')
ax4.set_ylabel(r'Nr. samples for SE $\leq 1$%')
ax2.set_ylabel(r'Nr. samples for SE $\leq 5$%')



h3, l3 = ax3.get_legend_handles_labels()
h4, l4 = ax4.get_legend_handles_labels()
ax3.legend(h3+h4,l3+l4)

plt.minorticks_off()
ax2.yaxis.set_minor_locator(MultipleLocator(1))
ax4.yaxis.set_major_locator(MultipleLocator(5))
ax4.yaxis.set_minor_locator(MultipleLocator(1))
ax3.yaxis.set_major_locator(MultipleLocator(0.01))


plt.tight_layout()
plt.savefig('size/size_required_samples.png',dpi=600)


fig,ax = plt.subplots(figsize=(4,3))

N_list = np.arange(1,101)


for i,size in enumerate([1.8, 2.3, 2.8, 3.5, 4.0, 5.0]):
    ax.plot(N_list,norm_error_comp111[i]/np.sqrt(N_list),label=f'{size} nm')

ax.set_xlabel('N samples')
ax.set_ylabel('111 composition SE')
ax.legend()

ax.set_ylim(0,0.15)

from matplotlib.ticker import MultipleLocator

ax.yaxis.set_major_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))

ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(5))

ax.grid(which='both',axis='y',alpha=0.8,ls=':')

plt.tight_layout()
plt.savefig('size/size_error_sample.png',dpi=600)



fig,ax = plt.subplots(figsize=(4,3))

N_list = np.arange(1,101)


for i,size in enumerate([1.8, 2.3, 2.8, 3.5, 4.0, 5.0]):
    ax.plot(N_list,df_std[i]/np.sqrt(N_list),label=f'{size} nm')

ax.set_xlabel('N samples')
ax.set_ylabel('$S_d$ SE')
ax.legend()

ax.set_ylim(0,0.05)

from matplotlib.ticker import MultipleLocator

ax.yaxis.set_major_locator(MultipleLocator(0.01))

ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(5))

ax.grid(which='both',axis='y',alpha=0.8,ls=':')

plt.tight_layout()
plt.savefig('size/size_error_sample_Sd.png',dpi=600)