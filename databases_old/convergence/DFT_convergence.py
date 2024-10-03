from ase.db import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


kpoints = np.arange(1,9)
with connect('kpoints_slab_out.db') as db, connect('kpoints_slab_def_out.db') as db_def:

    fcc111 = np.array([row.energy for row in db.select(surface='fcc111')])
    fcc100 = np.array([row.energy for row in db.select(surface='fcc100')])
    fcc211 = np.array([row.energy for row in db.select(surface='fcc211')])

    fcc111_def = np.array([row.energy for row in db_def.select(surface='fcc111')])
    fcc100_def = np.array([row.energy for row in db_def.select(surface='fcc100')])
    fcc211_def = np.array([row.energy for row in db_def.select(surface='fcc211')])



def plot_ax(ax,kpoints,E,surface):
    ax.scatter(kpoints,E)
    ax.set_xlabel('kpoints ([k,k,1])')
    ax.set_ylabel('$E_{slab} - E_{vacancy}$ [eV]')
    ax.set_title(surface)
    ax.set_xticks(kpoints)

    ax.axhline(np.mean(E[-3:]),ls='--',color='k',alpha=0.8,lw=1)


# fig,axes = plt.subplots(ncols=2,figsize=(7,6))

# for ax, E,surface in zip((ax1,ax2,ax3),(fcc111,fcc100,fcc211),('fcc111','fcc100','fcc211')):
#     ax.scatter(kpoints,E)
#     ax.set_xlabel('kpoints ([k,k,1])')
#     ax.set_ylabel('Total energy [eV]')
#     ax.set_title(surface)
#     ax.set_xticks(kpoints)


# plt.tight_layout()
# plt.savefig('kpoints_convergence.png',dpi=600)


# Plot 111 and 100
fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(6,3))
plot_ax(ax1,kpoints,fcc111-fcc111_def,'(111)')
plot_ax(ax2,kpoints,fcc100-fcc100_def,'(100)')

plt.tight_layout()
plt.savefig('kpoints_convergence1.png',dpi=600)
plt.close()

# Plot 211
fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(6,3))
plot_ax(ax1,kpoints,fcc211-fcc211_def,'(211)')

kpoints = np.arange(3,9)
with connect('kpoints_fcc211_out.db') as db,connect('kpoints_fcc211_def_out.db') as db_def:
    E = np.array([row.energy for row in db.select()]) - np.array([row.energy for row in db_def.select()])

plot_ax(ax2,kpoints,E,'(211)')
ax2.set_xlabel('kpoints ([k,4,1])')



# ax1.yaxis.set_major_locator(MultipleLocator(0.1))
# ax1.yaxis.set_minor_locator(MultipleLocator(0.025))

# ax2.yaxis.set_major_locator(MultipleLocator(0.05))
# ax2.yaxis.set_minor_locator(MultipleLocator(0.01))
# ax2.set_ylim(None,-248.548)

plt.tight_layout()
plt.savefig('kpoints_convergence2.png',dpi=600)
plt.close()



# fig,ax = plt.subplots()
# kpts=np.arange(3,9)
# ax4.scatter(kpts,E)
# ax4.set_xlabel('kpoints ([k,4,1])')
# ax4.set_ylabel('Total energy [eV]')
# ax4.set_title('fcc221')
# ax4.set_xticks(kpts)

# plt.tight_layout()
# plt.savefig('kpoints_convergence_fcc211_2.png',dpi=600)
# plt.tight_layout()
# plt.savefig('kpoints_convergence.png',dpi=600)



with connect('ecuts_slab_out.db') as db, connect('ecuts_slab_def_out.db') as db_def, connect('bulk_ecut.db') as db_bulk:

    fcc111 = np.array([row.energy for row in db.select(surface='fcc111')])
    fcc100 = np.array([row.energy for row in db.select(surface='fcc100')])
    fcc211 = np.array([row.energy for row in db.select(surface='fcc211')])

    fcc111_def = np.array([row.energy for row in db_def.select(surface='fcc111')])
    fcc100_def = np.array([row.energy for row in db_def.select(surface='fcc100')])
    fcc211_def = np.array([row.energy for row in db_def.select(surface='fcc211')])

    bulk = np.array([row.energy for row in db_bulk.select()])

    e = [row.ecut for row in db.select(surface='fcc111')]

print(e)
fig,axes = plt.subplots(ncols=3,figsize=(10,3))

for ax, E,surface in zip(axes,(fcc111-fcc111_def,fcc100-fcc100_def,fcc211-fcc211_def),('(111)','(100)','(211)')):
    ax.scatter(e,(E-bulk)*1000)
    ax.set_xlabel('$E_{cut}$ [meV]')
    ax.set_ylabel(r'$\Delta E$ [meV]')
    ax.set_title(surface)

    ax.set_xticks(ticks=e,labels=['200', '', '300', '', '400', '', '500', '600', '700', '800'])

    print([round((E[i]-E[-1])/E[-1],4) for i in range(len(E))])

plt.tight_layout()
plt.savefig('ecut_convergence.png',dpi=600)



