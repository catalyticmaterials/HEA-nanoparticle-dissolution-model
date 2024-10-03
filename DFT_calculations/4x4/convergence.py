from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import MultipleLocator

kpoints = np.arange(1,9)

fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(8,8))
with connect('Pt_kpoints_out.db') as db:
    
    for surface,ax in zip(('T111','T100','edge','kink'),axes.flatten()):
        slab = np.array([row.energy for row in db.select(surface=surface,defect='none')])
        slab_vac = np.array([row.energy for row in db.select(surface=surface,defect='vacancy')])

        ax.scatter(kpoints,slab_vac-slab)
        ax.set_title(surface)
        ax.set_xlabel('k-points [k,k,1]')

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        
        

plt.tight_layout()
plt.savefig('kpoints.png',dpi=600,bbox_inches='tight')
# plt.show()

# fig,ax = plt.subplots(figsize=(4,4))
# with connect('Pt_kpoints_out.db') as db:
    

#     slab = np.array([row.energy for row in db.select(surface='edge',defect='none')])
#     slab_vac = np.array([row.energy for row in db.select(surface='edge',defect='vacancy')])

#     ax.scatter(kpoints,slab_vac-slab)
#     ax.set_title('Edge')
#     ax.set_xlabel('k-points [k,4,1]')
#     ax.set_ylabel('$E_{vacancy} - E_{slab}$ [eV]')


# plt.tight_layout()
# plt.savefig('kpoints_edge.png',dpi=600,bbox_inches='tight')




fig,axes = plt.subplots(ncols=2,nrows=2,figsize=(8,8))
with connect('Pt_ecut_out.db') as db, connect('../Pt_ecut_bulks.db') as db_bulk:

    ecut = np.array([row.ecut for row in db_bulk.select()])
    bulk = np.array([row.energy for row in db_bulk.select()])
    
    for surface,ax in zip(('T111','T100','edge','kink'),axes.flatten()):
        slab = np.array([row.energy for row in db.select(surface=surface,defect='none')])
        slab_vac = np.array([row.energy for row in db.select(surface=surface,defect='vacancy')])
        
        print(slab,slab_vac)

        ax.scatter(ecut,slab_vac+bulk-slab)
        ax.set_title(surface)
        ax.set_xlabel('$E_{cut}$ [meV]')
        ax.set_ylabel(r'$\Delta E$ [eV]')
        

plt.tight_layout()
plt.savefig('ecut.png',dpi=600,bbox_inches='tight')