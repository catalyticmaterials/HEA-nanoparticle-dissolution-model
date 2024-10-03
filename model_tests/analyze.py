import numpy as np
import matplotlib.pyplot as plt
from utilities import metals

Ni = 1925


color_dict = {'batch_relax':'darkorange','batch_none':'tab:blue','single_relax':'forestgreen'}



fig, axes = plt.subplots(nrows=2,ncols=4,figsize=(12,6))
axes = axes.flatten()
for i,method in enumerate(('batch_none','batch_relax','single_relax')):
    
    data = np.loadtxt(f'{method}.csv',delimiter=',')

    
    surf_comp,surf111_comp,bulk_comp,surf_natoms, bulk_natoms, cn_dist,Sd,Nf = data[:,:8], data[:,8:16],data[:,24:32], data[:,32], data[:,33], data[:,34:44], data[:,44], data[:,45]

    sc_means = np.mean(surf_comp,axis=0)
    bc_means = np.mean(bulk_comp,axis=0)
    sc111_means = np.mean(surf111_comp,axis=0)
    sc_stds = np.std(surf_comp,axis=0,ddof=1)
    bc_stds = np.std(bulk_comp,axis=0,ddof=1)
    sc111_stds = np.std(surf111_comp,axis=0,ddof=1)

    diss_natoms = Ni - Nf    

    axes[0].errorbar(np.arange(8)-0.2+0.2*i,sc_means*100,sc_stds*100,alpha=1,fmt='.',color=color_dict[method],capsize=2)
    axes[0].set_xticks(range(8),metals)
    axes[0].set_title('Surface composition')


    axes[1].errorbar(np.arange(8)-0.2+0.2*i,sc111_means*100,sc111_stds*100,alpha=1,fmt='.',color=color_dict[method],capsize=2)
    axes[1].set_xticks(range(8),metals)
    axes[1].set_title('(111) surface composition')


    axes[2].errorbar(np.arange(8)-0.2+0.2*i,bc_means*100,bc_stds*100,alpha=1,fmt='.',color=color_dict[method],capsize=2)
    axes[2].set_xticks(range(8),metals)
    axes[2].set_title('Dissolved composition')
    
    cn_dist_norm = cn_dist/Nf.reshape(-1,1)
    cn_dist_norm[np.isnan(cn_dist_norm)] = 0.0
    axes[3].bar(np.arange(3,13)+i*0.2 - 0.3,np.mean(cn_dist_norm,axis=0) * 100,alpha=1,color=color_dict[method],label=method,width=0.2)
    axes[3].set_xticks(range(3,13))
    axes[3].set_title('Average number of each CN')

    if method!='batch_none':
        axes[4].hist(surf_natoms,bins=13,range=(400,465),alpha=0.6,color=color_dict[method])
        axes[4].set_title('$N_{atoms}$ in surface')
        axes[4].axvline(np.mean(surf_natoms),ls=':',color=color_dict[method])
        
        axes[6].hist(diss_natoms,bins=22,range=(360,580),alpha=0.6,color=color_dict[method])
        axes[6].set_title('$N_{atoms}$ dissolved')
        axes[6].axvline(np.mean(diss_natoms),ls=':',color=color_dict[method])

        axes[7].hist(Nf,bins=30,range=(1300,1600),alpha=0.6,color=color_dict[method])
        axes[7].set_title('$N_{atoms}$ total')
        axes[7].axvline(np.mean(Nf),ls=':',color=color_dict[method])

        axes[5].axvline(np.mean(Sd),ls=':',color=color_dict[method])
        axes[5].hist(Sd,bins=10,range=(0.25,0.45),alpha=0.6,color=color_dict[method])
        axes[5].set_title('(111) stability ($S_d$)')
    else:
        axins1 = axes[4].inset_axes((0.1,0.775,0.3,0.2))
        axins1.hist(surf_natoms,bins=30,range=None,alpha=0.6,color=color_dict[method])
        axins1.axvline(np.median(surf_natoms),ls=':',color=color_dict[method])
        axins1.minorticks_on()
        axins1.yaxis.set_tick_params(which='minor', bottom=False)

        axins2 = axes[6].inset_axes((0.65,0.775,0.3,0.2))
        axins2.hist(diss_natoms,bins=30,range=(900,2050),alpha=0.6,color=color_dict[method])
        axins2.axvline(np.median(diss_natoms),ls=':',color=color_dict[method])
        axins2.minorticks_on()
        axins2.yaxis.set_tick_params(which='minor', bottom=False)

        axins3 = axes[7].inset_axes((0.1,0.775,0.3,0.2))
        axins3.hist(Nf,bins=30,range=None,alpha=0.6,color=color_dict[method])
        axins3.axvline(np.median(Nf),ls=':',color=color_dict[method])
        axins3.minorticks_on()
        axins3.yaxis.set_tick_params(which='minor', bottom=False)

        axins4 = axes[5].inset_axes((0.65,0.775,0.3,0.2))
        axins4.hist(Sd,bins=15,range=(0,0.3),alpha=0.6,color=color_dict[method])
        axins4.axvline(np.median(Sd),ls=':',color=color_dict[method])
        axins4.minorticks_on()
        axins4.yaxis.set_tick_params(which='minor', bottom=False)

    

for ax in axes[4:]:
    ax.minorticks_on()
    ax.yaxis.set_tick_params(which='minor', bottom=False)



plt.tight_layout()

plt.subplots_adjust(bottom=0.1)
fig.legend(loc='outside lower center',ncols=3)


plt.savefig('test.png',dpi=600)

