import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
from utilities.compositionspace_functions import make_ternary_plot, get_molar_fractions, molar_fractions_to_cartesians
from utilities.stability import metals
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import itertools as it


# Get ternary grid
ternary_mfs = get_molar_fractions(0.1,3)
ternary_grid = molar_fractions_to_cartesians(ternary_mfs).T

# Set colormap
# sm = ScalarMappable(cmap='viridis')
sm = ScalarMappable(cmap='coolwarm_r')



for U in ('08','1'):

    data = np.loadtxt(f'grid_search_{U}V.csv',delimiter=',',skiprows=1)

    mf = data[:,:5]
    fd = data[:,5]

    # Make pseudo-ternary by combining IrRhRu
    # IrRhRu = np.sum(mf[:,[0,3,4]],axis=1)
    IrPt = np.sum(mf[:,[0,2]],axis=1).reshape(-1,1)
    RhRu = np.sum(mf[:,[3,4]],axis=1).reshape(-1,1)
    Pd = mf[:,1].reshape(-1,1)
    # p_ternary = np.hstack((IrRhRu.reshape(-1,1),mf[:,1:3]))
    p_ternary = np.hstack((IrPt,Pd,RhRu))
    # Get maximum of each duplicate
    p_fd = [np.max(fd[np.all(np.isclose(p_ternary,ternary_mf),axis=1)]) for ternary_mf in ternary_mfs]

    # Make pseudo_ternary_plot
    fig,ax=make_ternary_plot(ternary_grid,p_fd,['IrPt','Pd','RhRu'],vmin=0.0,vmax=1.0,colormap='coolwarm_r',minval=0.0)
    plt.colorbar(sm,ax=ax,shrink=0.5,anchor=(0.0,0.85),ticks=np.arange(0,1.1,0.2),label='$f_d$')

    plt.savefig(f'ternary_plots_{U}V/pseudo_ternary.png',dpi=600,bbox_inches='tight')
    plt.close()

    fig,axes = plt.subplots(nrows=2,ncols=5,figsize=(15,6))

    # Make each ternary plot
    idx_comb = it.combinations(np.arange(5),3)
    for ids,ax in zip(idx_comb,axes.flatten()):

        composition_mask = np.isclose(np.sum(mf[:,ids],axis=1),1)

        composition = [metals[i] for i in ids]

        alloy = ''.join(composition)

        ax = make_ternary_plot(ternary_grid,fd[composition_mask],composition,vmin=0.0,vmax=1.0,ax=ax,colormap='coolwarm_r',minval=0.0)
        
    
    plt.subplots_adjust(bottom=0.05)
    cbar_ax = fig.add_axes([0.2, 0.0, 0.6, 0.03])

    cbar = plt.colorbar(sm,cax=cbar_ax,ticks=np.arange(0,1.1,0.2),orientation='horizontal')
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('$f_d$',size=18)
    plt.tight_layout()
    plt.savefig(f'ternary_plots_{U}V/ternaries.png',dpi=600,bbox_inches='tight')
    plt.close()

