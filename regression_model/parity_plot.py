import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from math import log10,floor
from utilities.colors import metal_colors
from utilities import metals
from matplotlib.ticker import MultipleLocator


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        if num==0.0:
            return '0.0'
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)


def calc_R2(y_pred,y_true):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    # Calculate the total sum of squares (SS_tot)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    
    # Calculate the residual sum of squares (SS_res)
    ss_res = np.sum((y_true - y_pred) ** 2)
    
    # Calculate R-squared
    return  1 - (ss_res / ss_tot)
    



def partity_plot(y_true,y_preds,metal_feature,unit='eV'):
    y_true = np.array(y_true)
    y_preds = np.array(y_preds)
    properties = {'eV':r'$\Delta E$','V':'U'}
    prop = properties[unit]

    E = y_preds - y_true
    MAE = np.mean(np.abs(E))
    RMSE = np.sqrt(np.mean(E**2))
    R2 = calc_R2(y_preds,y_true)
    print('RMSE:',RMSE)
    print('R2:',R2)


    fig, ax = plt.subplots(figsize=(3.,3.))

    # Get MAE of each metal
    print('metal: MAE, MASE')
    labels = {}
    for i,metal in enumerate(metals):
        mask = metal_feature[:,i]==1
        MAE_metal = np.mean(np.abs(E[mask]))
        MAD_metal = np.mean(np.abs(y_true[mask]-np.mean(y_true[mask])))
        MASE = MAE_metal/MAD_metal


        labels[metal] = f'{metal}: {MAE_metal:1.3f}({MASE:1.2f})'
        print(metal,MAE_metal,MASE)


    ax_ins = ax.inset_axes((0.51,0.18,0.48,0.38))


    colors = [metal_colors[metals[list(x).index(1)]] for x in metal_feature]


    ax.scatter(y_true,y_preds,marker='.',color=colors,zorder=1,alpha=0.7,edgecolors='none',s=50)

    lims = np.array([np.min([y_true,y_preds])-0.2,np.max([y_true,y_preds])+0.2])


    ax.plot(lims,lims,c='k',zorder=0)


    ax_ins.hist(E,bins=int(np.sqrt(len(E))),color='k',alpha=0.33)
    for i,metal in enumerate(metals):
        ax_ins.hist(E[metal_feature[:,i]==1],bins=15,color=metal_colors[metal],alpha=0.7,histtype='step')
    ax_ins.axvline(np.mean(E),c='k',ls=':',lw=1)
    # Remove borders, labels, and ticks
    ax_ins.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax_ins.spines['right'].set_visible(False)
    ax_ins.spines['top'].set_visible(False)
    ax_ins.spines['left'].set_visible(False)
    ax_ins.yaxis.tick_left()
    ax_ins.tick_params(axis='y', left=False, labelleft=False)
    ax_ins.tick_params(axis='x', labelsize=8)
    ax_ins.set_xlabel(f'{prop}$_{{pred}} -$ {prop}$_{{DFT}}$' + f'[{unit}]',fontsize=8)
    ax_ins.set_facecolor('none')

    ax.set_ylim(lims)
    ax.set_xlim(lims)


    ax.set_xlabel(f'{prop}$_{{DFT}}$ [{unit}]')
    ax.set_ylabel(f'{prop}$_{{pred}}$ [{unit}]')


    txt = f' All: {MAE:1.3f}'
    ax.text(0.0055,0.54,txt,va='top',transform=ax.transAxes,fontsize=8,family='Consolas')
    
    handles = []
    for metal in metals:
        handles.append(Line2D([0], [0], marker='.', color="w", label=labels[metal],markerfacecolor=metal_colors[metal])) 


    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    

    L = ax.legend(handles=handles,fontsize=8,title_fontproperties={'family': 'Consolas','size':8},loc='upper left',bbox_to_anchor=(0.01,0.995),alignment='left',borderpad=0.0,frameon=False,title=f'Metal: MAE [{unit}] (MASE)',labelspacing=0.1,markerscale=2.25,borderaxespad=0.0,handletextpad=0.0,handlelength=1,framealpha=1)
    plt.setp(L.texts, family='Consolas')
    ax.set_aspect('equal')
    return fig, ax