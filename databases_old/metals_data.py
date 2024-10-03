from ase.db import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from utilities.colors import metal_colors
from utilities.stability import metals, U_standard
from matplotlib.ticker import MultipleLocator

markers = ['s','o','^','X','d']



fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(6.5,3),sharex=True)
data_U = []
data_dE = []

with connect('metals_dissolution_out.db') as db, connect('bulks.db') as db_bulk:

    E_bulk = {row.metal:row.energy for row in db_bulk.select()}


    for metal,marker in zip(metals,markers):

        no_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='none')}

        add_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='add')}

        remove_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='removed')}


        # Delta energies of coordination numbers
        # dE_cn3 = add_defect['T111'] + E_bulk[metal]  - no_defect['T111']

        # dE_cn4 = add_defect['T100'] + E_bulk[metal]  - no_defect['T100']

        # dE_cn5 = add_defect['edge'] + E_bulk[metal]  - no_defect['edge']

        # dE_cn6 = remove_defect['kink'] + E_bulk[metal] -no_defect['kink']

        # dE_cn7 = remove_defect['edge'] + E_bulk[metal] -no_defect['edge']

        # dE_cn8 = remove_defect['T100'] + E_bulk[metal] -no_defect['T100']

        # dE_cn9 = remove_defect['T111'] + E_bulk[metal] -no_defect['T111']


        dE_add = [no_defect[surface] + E_bulk[metal]  - add_defect[surface] for surface in ('T111','T100','edge')]

        dE_remove = [remove_defect[surface] + E_bulk[metal] -no_defect[surface] for surface in ('kink','edge','T100','T111')]

        dE = dE_add + dE_remove
        dE = np.array(dE)

        ax1.scatter(range(3,10),dE, color=metal_colors[metal],marker=marker)

        U_metal = U_standard[metal]
        
        U = dE/U_metal['n'] + U_metal['U']

        ax2.scatter(range(3,10),U, color=metal_colors[metal],marker=marker)

        data_dE.append(dE)
        data_U.append(U)


data_dE = np.array(data_dE).T
data_U = np.array(data_U).T

data = np.hstack((data_dE,data_U))
np.savetxt('metals_data.csv',data,delimiter=',',header='dE_Ir,dE_Pd,dE_Pt,dE_Rh,dE_Ru,U_Ir,U_Pd,U_Pt,U_Rh,U_Ru',fmt='%1.6f')



handles = []
for metal,marker in zip(metals,markers):
    handles.append(Line2D([0], [0], marker=marker, color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=12)) 

ax1.set_ylabel('$\Delta E$ [eV]')
ax2.set_ylabel('Dissolution potential [V]')



for ax in (ax1,ax2):
    ax.set_xlabel('Coordination Number')
    ax.set_xticks(np.arange(3,10))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))


pos1 = ax1.get_position()
pos2 = ax2.get_position()

fig.legend(handles=handles, labels=metals,
           loc='outside upper center', ncol=5, mode='expand',fontsize=12,bbox_to_anchor=(pos1.x0, .5, pos2.x1-pos1.x0, 0.5),fancybox=False)


plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.savefig('Metals.png',dpi=600,bbox_inches='tight')
