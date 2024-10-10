import numpy as np
import matplotlib.pyplot as plt
from utilities.colors import alloy_color
from utilities.compositionspace_functions import prepare_triangle_plot, get_molar_fractions, molar_fractions_to_cartesians, make_ternary_plot
from matplotlib.cm import ScalarMappable


mfs = get_molar_fractions(0.01,3)
xgrid = molar_fractions_to_cartesians(mfs)


for ternary_metals in (['Pd','Pt','Ru'],['Au','Cu','Pt'],['Au','Cu','Pd']):
    break
    system = ''.join(ternary_metals)

    c = [alloy_color(ternary_metals,mf) for mf in mfs]

    fig,ax = plt.subplots(figsize=(4,4))
    ax = prepare_triangle_plot(ax,ternary_metals,edges=False,ticks=False)
    ax.scatter(xgrid.T[0],xgrid.T[1],color=c,marker='h',s=1.85)

    plt.savefig(f'ternaries/{system}_cmap.png',dpi=600,bbox_inches='tight')

    data=np.loadtxt(f'ternaries/{system}_sc_maping.csv',delimiter=',')

    bulk_comp = data[:,:3]
    isc = data[:,3:6]
    fsc = data[:,6:-1]
    Sd = data[:,-1]

    bulk_colors = np.array([alloy_color(ternary_metals,mf) for mf in bulk_comp])

    # stab_mask = np.isclose(np.sum(fsc,axis=1),1)
    stab_mask = Sd>0.0
    unstab_mask = np.invert(stab_mask)

    fsc[unstab_mask] = 0.0

    fig,ax = plt.subplots(figsize=(4,4))
    ax = prepare_triangle_plot(ax,ternary_metals)
    
    X = molar_fractions_to_cartesians(fsc[stab_mask]).T
    X_unstab = molar_fractions_to_cartesians(bulk_comp[unstab_mask]).T

    for i,x in enumerate(X.T):
        ax.scatter(x[0],x[1],color=bulk_colors[stab_mask][i],marker='o',s=50*Sd[stab_mask][i]+10,alpha=0.8,edgecolors='none')
    ax.scatter(X_unstab[0],X_unstab[1],facecolors='none',edgecolors=bulk_colors[unstab_mask],marker='o')

    X_isc = molar_fractions_to_cartesians(isc[stab_mask])
    for i in range(X.shape[1]):
        ax.plot([X_isc[i,0],X[0,i]],[X_isc[i,1],X[1,i]],c='k',alpha=0.08,zorder=0)


    plt.savefig(f'ternaries/{system}_fsc2bc.png',dpi=600,bbox_inches='tight')


    fsc_colors = [alloy_color(ternary_metals,mf) for mf in fsc]

    fig,ax = plt.subplots(figsize=(4,4))
    ax = prepare_triangle_plot(ax,ternary_metals)
    
    X = molar_fractions_to_cartesians(bulk_comp).T

    ax.scatter(X[0],X[1],color=fsc_colors,marker='h',s=115)

    plt.savefig(f'ternaries/{system}_bc2fsc.png',dpi=600,bbox_inches='tight')



    fig,ax = plt.subplots(figsize=(4,4))
    ax=make_ternary_plot(X,Sd,ternary_metals,ax=ax,vmin=0.0,vmax=1.0,colormap='coolwarm_r',minval=0.0)

    plt.savefig(f'ternaries/{system}_Sd.png',dpi=600,bbox_inches='tight')