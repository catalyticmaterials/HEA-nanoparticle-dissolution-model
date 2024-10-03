from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable


kpoints = np.arange(1,9)

fig,axes = plt.subplots(2,2,figsize=(6,5),sharex=True)
# plt.axis('off')
# for ax in axes.flatten():
#     ax.spines[['left','bottom']].set_visible(False)
#     ax.set_xticks([],[])
#     ax.set_yticks([],[])
# plt.subplots_adjust(hspace=0.45,wspace=0.4)


# sps1,sps2,sps3,sps4 = GridSpec(2,2,fig)

# bax1 = brokenaxes(ylims=((-1.04,-1.02),(-0.015,0.02)), subplot_spec=sps1,hspace=0.2)
# bax2 = brokenaxes(ylims=((-0.7,-0.5),(-0.2,0.1)), subplot_spec=sps2,hspace=0.2)
# bax3 = brokenaxes(ylims=((-0.32,-0.25),(-0.05,0.1)), subplot_spec=sps3,hspace=0.2)
# bax4 = brokenaxes(subplot_spec=sps4)
# bax1.set_yticks([-1.04,-1.02])

# baxes = [bax1,bax2,bax3,bax4]
# baxes=[bax1,axes[0,1],axes[1,0],axes[1,1]]

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

with connect('Pt_kpoints_out.db') as db:
    
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
        


        # d = .5  # proportion of vertical to horizontal extent of the slanted line
        # kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
        #             linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        # ax.plot([0, 1], [0, 0], transform=ax.transAxes, **kwargs)
        # ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
        


        
        # ax.set_xlabel('k-points [k,k,1]')
        # ax.set_ylabel('$E_{vacancy} - E_{slab}$ [eV]')
        # ax.set_ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]')

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(locator[0]))
        ax.yaxis.set_minor_locator(MultipleLocator(locator[1]))
        
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)
        # break
        # ax.spines['right'][0].set_visible(False)
# print(bax1.spines)
# axes[0,0].yaxis.set_major_locator(MultipleLocator(0.5))

# ylim=axes[0,0].get_ylim()
# breaky = ((ylim[0],-0.8),(-0.2,ylim[1]))
# print(breaky)
# bax1.spines['right'][0].set_visible(True)
# bax1.spines['top'][0].set_visible(True)
# axes[0,0].tick_params('both','both',)

# axes[0,0].xaxis.set_visible(False)
# axes[0,0].yaxis.set_visible(False)
# axes[0,0].set_ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]')
# axes[0,0].spines['left'].set_visible(False)
# bax1.set_ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]')

axes[1,0].set_xlabel('k-points [k,k,1]')
axes[1,1].set_xlabel('k-points [k,k,1]')
# axes[0,0].set_ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]',loc='top',labelpad=10)
# axes[1,0].set_ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]')

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# plt.xlabel("common X")
plt.ylabel(r'$\Delta E(k) - \Delta E(k=8)$ [eV]',labelpad=20)

# plt.tight_layout()
plt.subplots_adjust(wspace=0.4,hspace=0.25)

plt.savefig('kpoints.png',dpi=600,bbox_inches='tight')
# plt.show()

fig,ax = plt.subplots(figsize=(4,4))
with connect('Pt_kpoints_out.db') as db:
    

    slab = np.array([row.energy for row in db.select(surface='edge',defect='none')])
    slab_vac = np.array([row.energy for row in db.select(surface='edge',defect='vacancy')])

    ax.scatter(kpoints,slab_vac-slab)
    ax.set_title('Edge')
    ax.set_xlabel('k-points [k,4,1]')
    ax.set_ylabel('$E_{vacancy} - E_{slab}$ [eV]')


plt.tight_layout()
plt.savefig('kpoints_edge.png',dpi=600,bbox_inches='tight')




fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(6,5),sharex=True)

with connect('Pt_ecut_out.db') as db, connect('Pt_ecut_bulks.db') as db_bulk:

    ecut = np.array([row.ecut for row in db_bulk.select()])
    bulk = np.array([row.energy for row in db_bulk.select()])
    
    for surface,ax in zip(('T111','T100','edge','kink'),axes.flatten()):
        slab = np.array([row.energy for row in db.select(surface=surface,defect='none')])
        slab_vac = np.array([row.energy for row in db.select(surface=surface,defect='vacancy')])

        dE = slab_vac+bulk-slab
        dE -=dE[-1]
        dE*=1000
        # if surface=='kink':
            # dE*=1000
            # ax.set_ylabel(r'$\Delta E(E_{cut}) - \Delta E(E_{cut}=800)$ [meV]')
        #     pass
        # else:
            # ax.set_ylabel(r'$\Delta E(E_{cut}) - \Delta E(E_{cut}=800)$ [eV]')
            # ax.yaxis.set_major_locator(MultipleLocator(0.01))
        # ax.set_ylabel(r'$\Delta E(E_{cut}) - \Delta E(E_{cut}=800)$ [eV]')

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
# plt.xlabel("common X")
plt.ylabel(r'$\Delta E(E_{cut}) - \Delta E(E_{cut}=800)$ [meV]',labelpad=15)


plt.tight_layout()
plt.savefig('ecut.png',dpi=600,bbox_inches='tight')