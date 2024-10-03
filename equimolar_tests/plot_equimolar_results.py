import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from utilities import metals
from utilities.colors import metal_colors
from scipy.stats import norm,kstest


def kstest_normal(dist):
    mean,std = np.mean(dist), np.std(dist,ddof=1)
    pnorm = kstest(dist,norm(mean,std).cdf).pvalue
    return mean,std,pnorm

n_metals = len(metals)


# rcParams.update({'font.size': 12})
for size in ('35A','4nm'):
    data = np.loadtxt(f'equimolar{size}_muM_results.csv',delimiter=',',skiprows=1)
    final_111_compositions = data[:,:n_metals]
    dissolution_compositions = data[:,n_metals:-1]
    dissolution_factors = data[:,-1]


    bins = np.arange(-0.005,np.round(np.max(final_111_compositions),decimals=2)+0.005,0.01)
    fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4.2),sharey=True,sharex=True)
    axes[0,0].set_yticks([0,100,200],[0,50,100])
    axes[0,0].set_ylabel('Frequency [%]')
    axes[1,0].set_ylabel('Frequency [%]')
    axes[0,0].set_ylim(0.0,200)
    axes[0,0].yaxis.set_minor_locator(MultipleLocator(20))
    axes[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    axes = axes.flatten()
    for i, metal in enumerate(metals):

        surface = final_111_compositions[:,i]

        axes[i].hist(surface,color=metal_colors[metal],bins=bins)
  
        axes[i].set_xlabel(f'{metal}$_x$ Surface')

       
        mean,std,pnorm = kstest_normal(surface)
        text = f'$\mu={mean:1.2f}, \sigma={std:1.2f}$\n$p_{{norm}}={pnorm:1.2f}$'
        axes[i].text(0.98,0.98,text,ha='right',va='top',transform=axes[i].transAxes)


    


    plt.tight_layout()
    plt.savefig(f'Equimolar{size}_muM_dissolution_111_compositions.png',dpi=600,bbox_inches='tight')
    plt.close()


    fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharey=True,sharex=True)
    axes[0,0].set_yticks([0,100,200],[0,50,100])
    axes[0,0].set_ylabel('Frequency [%]')
    axes[1,0].set_ylabel('Frequency [%]')
    axes[0,0].set_ylim(0.0,200)
    axes[0,0].yaxis.set_minor_locator(MultipleLocator(20))
    axes[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    axes = axes.flatten()
    bins = np.arange(-0.005,np.round(np.max(dissolution_compositions),decimals=2)+0.005,0.01)
    for i, metal in enumerate(metals):

        diss = dissolution_compositions[:,i]

        axes[i].hist(diss,color=metal_colors[metal],bins=bins)
        
        axes[i].set_xlabel(f'{metal}$_x$ Dissolved')

        mean,std,pnorm = kstest_normal(diss)
        if np.isnan(pnorm):
            pnorm=0.0
        text = f'$\mu={mean:1.2f}, \sigma={std:1.2f}$\n$p_{{norm}}={pnorm:1.2f}$'
        axes[i].text(0.98,0.98,text,ha='right',va='top',transform=axes[i].transAxes)
        



    plt.tight_layout()
    plt.savefig(f'Equimolar{size}_muM_dissolution_diss_compositions.png',dpi=600,bbox_inches='tight')
    plt.close()




    fig,ax = plt.subplots(figsize=(4,3))

    ax.hist(dissolution_factors,color='grey',bins=12)
    ax.set_xlabel('$f_d$')
    ax.set_ylabel('Frequency')

    mean,std,pnorm = kstest_normal(dissolution_factors)
    text = f'$\mu={mean:1.2f}, \sigma={std:1.2f}$\n$p_{{norm}}={pnorm:1.2f}$'
    ax.text(0.98,0.98,text,ha='right',va='top',transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(f'Equimolar{size}_muM_dissolution_factors.png',dpi=600,bbox_inches='tight')
    break