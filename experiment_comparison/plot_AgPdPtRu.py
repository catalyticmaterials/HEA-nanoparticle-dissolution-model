import matplotlib.pyplot as plt
import numpy as np
from utilities. compositionspace_functions import prepare_triangle_plot, molar_fractions_to_cartesians
from utilities.colors import alloy_color


all_metals = ['Ag','Pd','Pt','Ru']
data = np.loadtxt('AgPdPtRu1M.csv',delimiter=',',skiprows=1)
bulk_comp = data[:,:4]
surface_comp = data[:,4:8]
Sd = data[:,-1]

for metals in (['Ag','Pd','Ru'],['Ag','Pt','Ru']):

    ternary_mask = np.array([metal in metals for metal in all_metals])
    bulk_ternary = bulk_comp[:,ternary_mask]
    surface_ternary = surface_comp[:,ternary_mask]
    alloy_mask = np.isclose(np.sum(bulk_ternary,axis=1),1.0)

    
    bulk_mf = bulk_ternary[alloy_mask]
    surf_mf = surface_ternary[alloy_mask]

    surf_x = molar_fractions_to_cartesians(surf_mf)
    bulk_x = molar_fractions_to_cartesians(bulk_mf)
    bulk_color = np.array([alloy_color(metals,bulk_mf_i) for bulk_mf_i in bulk_mf])
    
    # non_zero_mask = np.invert(np.all(np.isclose(surf_mf,0.0),axis=1))
    non_zero_mask = Sd[alloy_mask]>0.0
    zero_mask = np.invert(non_zero_mask)

    fig,ax = plt.subplots(figsize=(4,4))
    ax = prepare_triangle_plot(ax,metals)

    ax.scatter(*surf_x[non_zero_mask].T,color=bulk_color[non_zero_mask])
    ax.scatter(*bulk_x[zero_mask].T,color=bulk_color[zero_mask],facecolor='none')
    for i in range(len(surf_x)):
        ax.plot([bulk_x[i,0],surf_x[i,0]],[bulk_x[i,1],surf_x[i,1]],c='k',alpha=0.25)

    alloy = ''.join(metals)
    plt.savefig(f'{alloy}_bc2fsc.png',dpi=600,bbox_inches='tight')


  
