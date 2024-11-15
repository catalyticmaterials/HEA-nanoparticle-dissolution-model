from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable


# k-points
kpoints = np.arange(1,9)

fig,axes = plt.subplots(2,2,figsize=(6,5),sharex=True)

titles = {'T111':'Terrace (111)','T100':'Terrace (100)','edge':'Edge','kink':'Kink'}

def break_y(ax,x,y,lim1,lim2):
    divider = make_axes_locatable(ax)
    ax2 = divider.new_vertical(size="300%", pad=0.25)
    fig.add_axes(ax2)

    ax.scatter(x,y)
    ax2.scatter(x,y)

    ax.set_ylim(*lim1)
    ax2.set_ylim(*lim2)

    ax.spines['top'].set_visible(False)
    ax2.tick_params(bottom=False, labelbottom=False)
    ax2.spines['bottom'].set_visible(False)

    d = .025  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False,linewidth=1)
    ax2.plot((-d, +d), (-d/3, +d/3), **kwargs)        # top-left diagonal
    ax2.plot((1 - d, 1 + d), (-d/3, +d/3), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
    ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    return ax,ax2

ylims = [((-1.04,-1.03),(-0.01,0.01)),((-0.65,-0.6),(-0.05,0.1)),((-0.32,-0.25),(-0.05,0.05)),None]

locators = [(0.01,0.001),(0.05,0.001),(0.05,0.001),(0.01,0.001)]

with connect('convergence/Pt_kpoints_out.db') as db:
    
    for surface,ax,lims,locator in zip(('T111','T100','edge','kink'),axes.flatten(),ylims,locators):
        slab = np.array([row.energy for row in db.select(surface=surface,defect='none')])
        slab_vac = np.array([row.energy for row in db.select(surface=surface,defect='vacancy')])

        dE = slab_vac-slab
        dE = dE - dE[-1]

        if surface=='kink':
            ax.scatter(kpoints,dE)
            ax.set_title(titles[surface])
        else:
            ax,ax2 = break_y(ax,kpoints,dE,lims[0],lims[1])
            ax2.yaxis.set_major_locator(MultipleLocator(locator[0]))
            ax2.yaxis.set_minor_locator(MultipleLocator(locator[1]))
            ax2.minorticks_on()
            ax2.xaxis.set_tick_params(which='minor', bottom=False)
            ax2.set_title(titles[surface])
    

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(locator[0]))
        ax.yaxis.set_minor_locator(MultipleLocator(locator[1]))
        
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)


axes[1,0].set_xlabel('k-points [k,k,1]')
axes[1,1].set_xlabel('k-points [k,k,1]')


# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# plt.xlabel("common X")
plt.ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]',labelpad=20)

# plt.tight_layout()
plt.subplots_adjust(wspace=0.4,hspace=0.25)

plt.savefig('convergence/kpoints.png',dpi=600,bbox_inches='tight')



# Ecut
fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(6,5),sharex=True)

with connect('convergence/Pt_ecut_out.db') as db, connect('convergence/Pt_ecut_bulks.db') as db_bulk:

    ecut = np.array([row.ecut for row in db_bulk.select()])
    bulk = np.array([row.energy for row in db_bulk.select()])
    
    for surface,ax in zip(('T111','T100','edge','kink'),axes.flatten()):
        slab = np.array([row.energy for row in db.select(surface=surface,defect='none')])
        slab_vac = np.array([row.energy for row in db.select(surface=surface,defect='vacancy')])

        dE = slab_vac+bulk-slab
        dE -=dE[-1]
        dE*=1000


        ax.scatter(ecut,dE)
        ax.set_title(titles[surface])
        
        
        ax.set_xticks([200,250,300,350,400,450,500,600,700,800],['200','','300','','400','','500','600','700','800'])

        
        # ax.yaxis.set_minor_locator(MultipleLocator(0.001))
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)

axes[1,0].set_xlabel('$E_{cut}$ [eV]') 
axes[1,1].set_xlabel('$E_{cut}$ [eV]') 

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

plt.ylabel(r'$\Delta E(E_{cut}) - \Delta E(E_{cut}=800)$ [meV]',labelpad=15)


plt.tight_layout()
plt.savefig('convergence/ecut.png',dpi=600,bbox_inches='tight')