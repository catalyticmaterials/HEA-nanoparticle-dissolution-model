from ase.db import connect
import numpy as np
import matplotlib.pyplot as plt
from utilities.colors import metal_colors
from utilities import metals, U_standard
from matplotlib.ticker import MultipleLocator

markers = ['P','p','H','s','o','^','X','d']



c=1e-6
kB=8.617333262e-5
T=298.15
U_dict = U_standard
for metal in metals:
    U_dict[metal]['U'] = U_dict[metal]['U'] + kB*T/U_dict[metal]['n'] * np.log(c)

# Extract data
data_U = []
data_dE = []
with connect('metal_dbs/metals_dissolution_out.db') as db, connect('metal_dbs/bulks.db') as db_bulk:

    E_bulk = {row.metal:row.energy for row in db_bulk.select()}


    for metal,marker in zip(metals,markers):

        no_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='none')}

        add_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='adatom')}

        remove_defect = {row.surface:row.energy for row in db.select(metal=metal,defect='vacancy')}



        dE_add = [no_defect[surface] + E_bulk[metal]  - add_defect[surface] for surface in ('T111','T100','edge')]

        dE_remove = [remove_defect[surface] + E_bulk[metal] -no_defect[surface] for surface in ('kink','edge','T100','T111')]

        dE = dE_add + dE_remove
        dE = np.array(dE)


        U_metal = U_dict[metal]
        
        U = np.min(dE/U_metal['n'] + U_metal['U'],axis=0)


        data_dE.append(dE)
        data_U.append(U)


data_dE = np.array(data_dE).T
data_U = np.array(data_U).T

data = np.hstack((data_dE,data_U))
np.savetxt('metals_data.csv',data,delimiter=',',header='dE_Ag,dE_Au,dE_Cu,dE_Ir,dE_Pd,dE_Pt,dE_Rh,dE_Ru,U_Ag,U_Au,U_Cu,U_Ir,U_Pd,U_Pt,U_Rh,U_Ru',fmt='%1.6f')




# Plot Delta E
fig,axes = plt.subplots(nrows=2,ncols=4,figsize=(8,4),sharex=True,sharey=True)
axes = axes.flatten()
for i,metal in enumerate(metals):

    axes[i].scatter(range(3,10),data_dE.T[i],color=metal_colors[metal],label=metal)

    axes[i].legend(loc='lower right')

    axes[i].set_xticks(range(3,10))

    if i>3:
        axes[i].set_xlabel('CN')

axes[0].set_ylabel(r'$\Delta E$ [eV]')
axes[4].set_ylabel(r'$\Delta E$ [eV]')
axes[0].yaxis.set_minor_locator(MultipleLocator(0.1))
axes[4].yaxis.set_minor_locator(MultipleLocator(0.1))

plt.tight_layout()
plt.savefig('metals_dE.png',dpi=600,bbox_inches='tight')




#  Plot Udiss
fig,ax = plt.subplots(figsize=(5,3))
for i,metal in enumerate(metals):

    ax.scatter(range(3,10),data_U.T[i],color=metal_colors[metal],marker=markers[i],label=metal)



ax.set_xticks(range(3,10))
ax.set_xlabel('CN')
ax.set_ylabel('$U_{dissolution}$ [V]')
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

plt.legend(loc='center left',bbox_to_anchor=(1.02,0.5),fancybox=False)


plt.tight_layout()
plt.savefig('metals_Udiss.png',dpi=600,bbox_inches='tight')

