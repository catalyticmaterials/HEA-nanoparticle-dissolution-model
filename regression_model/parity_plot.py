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

    properties = {'eV':r'$\Delta E$','V':'U'}
    prop = properties[unit]

    E = y_preds - y_true
    MAE = np.mean(np.abs(E))
    ME = np.mean(E)
    ME_SE = np.std(E,ddof=1)/np.sqrt(len(E))
    RMSE = np.sqrt(np.mean(E**2))
    # std = np.std(E)
    R2 = calc_R2(y_preds,y_true)

    # Get MAE of each metal
    MAE_metals_text = [f'Metal: MAE [{unit}]']
    for i,metal in enumerate(metals):
        mask = metal_feature[:,i]==1
        MAE_metal = np.mean(np.abs(E[mask]))
        MAE_metals_text.append(f'{metal}: {MAE_metal:1.3f}')

    MAE_metals_text = '\n'.join(MAE_metals_text)

    fig, ax = plt.subplots(figsize=(4.5,4.5))

    ax_ins = ax.inset_axes((0.51,0.15,0.48,0.38))


    colors = [metal_colors[metals[list(x).index(1)]] for x in metal_feature]


    ax.scatter(y_true,y_preds,marker='.',color=colors,zorder=1,alpha=0.7,edgecolors='none',s=50)

    lims = np.array([np.min([y_true,y_preds])-0.1,np.max([y_true,y_preds])+0.1])


    ax.plot(lims,lims,c='k',zorder=0)
    ax.plot(lims,lims - 0.1,c='k',ls='--',zorder=0)
    ax.plot(lims,lims + 0.1,c='k',ls='--',zorder=0)


    ax_ins.hist(E,bins=int(np.sqrt(len(E))),color='k',alpha=0.33)
    for i,metal in enumerate(metals):
        ax_ins.hist(E[metal_feature[:,i]==1],bins=15,color=metal_colors[metal],alpha=0.7,histtype='step')
    ax_ins.axvline(np.mean(E),c='k',ls=':',lw=1)
    # Remove borders, labels, and ticks
    ax_ins.spines['right'].set_visible(False)
    ax_ins.spines['top'].set_visible(False)
    ax_ins.spines['left'].set_visible(False)
    ax_ins.yaxis.tick_left()
    ax_ins.tick_params(axis='y', left=False, labelleft=False)
    ax_ins.set_xlabel(f'{prop}$_{{pred}} -$ {prop}$_{{DFT}}$' + f'[{unit}]')
    ax_ins.set_facecolor('none')

    ax.set_ylim(lims)
    ax.set_xlim(lims)


    ax.set_xlabel(f'{prop}$_{{DFT}}$ [{unit}]')
    ax.set_ylabel(f'{prop}$_{{pred}}$ [{unit}]')


    txt = f'MAE = {MAE:1.3f} {unit}\nRMSE = {RMSE:1.3f} {unit}\n$R^2$ = {R2:1.3f} {unit}'
    ax.text(0.23,0.975,txt,va='top',transform=ax.transAxes)
    
    ax.text(0.23,0.825,MAE_metals_text,va='top',transform=ax.transAxes,fontsize=8)

    handles = []
    for metal in metals:
        handles.append(Line2D([0], [0], marker='o', color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=10)) 


    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax_ins.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.legend(handles=handles,loc=2)
    ax.set_aspect('equal')
    return fig, ax