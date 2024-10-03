from ase.db import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from utilities.stability import metals, U_standard
from utilities.colors import metal_colors



with connect('bulks.db') as db:
    E_bulk={row.metal:row.energy for row in db.select()}


# Number of different slab compositions
N = 50



# Function to handle 'no match' error in database
def util_func(db,slab_idx,metal):
    try: 
        return db.get(slab_idx=slab_idx,metal=metal,defect='replace').energy
    except KeyError:
        return db.get(slab_idx=slab_idx,metal=metal,defect='none').energy


def ax_violin(ax1,ax2,dE,U,cn,metal):
    parts1 = ax1.violinplot(dE,[cn],showextrema=False)
    parts2 = ax2.violinplot(U,[cn],showextrema=False)

    for pc in parts1['bodies']:
        pc.set_facecolor(metal_colors[metal])
        pc.set_edgecolor(metal_colors[metal])
        pc.set_alpha(1)

    for pc in parts2['bodies']:
        pc.set_facecolor(metal_colors[metal])
        pc.set_edgecolor(metal_colors[metal])
        pc.set_alpha(1)



def get_dE_and_U(db_path,adatom=True):
    with connect(db_path) as db: 
        
        dE_rem = np.array([db.get(slab_idx=i,metal='none',defect='removed').energy + E_bulk[metal] - util_func(db,i,metal) for i in range(N)])
        if adatom:
            dE_ad = np.array([db.get(slab_idx=i,defect='none').energy + E_bulk[metal] - db.get(slab_idx=i,metal=metal,defect='add').energy for i in range(N)])

    U_rem = dE_rem/U_metal['n'] + U_metal['U']

    if adatom:
        U_ad = dE_ad/U_metal['n'] + U_metal['U']
        return dE_rem, U_rem, dE_ad, U_ad
    
    return dE_rem, U_rem
    
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

fig, axes = plt.subplots(ncols=5,nrows=2,figsize=(15,6),sharex=True)
# ((ax11,ax12,ax13,ax14,ax15),(ax21,ax22,ax23,ax24,ax25))
for metal, (ax1,ax2) in zip(metals,list(zip(*axes))):

    U_metal = U_standard[metal]


    # Get T111 data
    dE_111_rem,U_111_rem,dE_111_add,U_111_add = get_dE_and_U('hea_dbs/T111_HEA_out.db')

    ax_violin(ax1,ax2,dE_111_add,U_111_add,3,metal)
    ax_violin(ax1,ax2,dE_111_rem,U_111_rem,9,metal)
    

    # Get T100 data
    dE_100_rem,U_100_rem, dE_100_add,U_100_add = get_dE_and_U('hea_dbs/T100_HEA_out.db')

    ax_violin(ax1,ax2,dE_100_add,U_100_add,4,metal)
    ax_violin(ax1,ax2,dE_100_rem,U_100_rem,8,metal)


    # Get edge data
    dE_edge_rem,U_edge_rem, dE_edge_add,U_edge_add = get_dE_and_U('hea_dbs/Edge_HEA_out.db')

    ax_violin(ax1,ax2,dE_edge_add,U_edge_add,5,metal)
    ax_violin(ax1,ax2,dE_edge_rem,U_edge_rem,7,metal)


    # Get kink data
    dE_kink_rem,U_kink_rem = get_dE_and_U('hea_dbs/Kink_HEA_out.db',adatom=False)

    ax_violin(ax1,ax2,dE_kink_rem,U_kink_rem,6,metal)



    # ax1.set_xlabel('Coordination Number')
    # ax2.set_xlabel('Coordination Number')
    # ax1.set_ylabel('$\Delta E$ [eV]')
    # ax2.set_ylabel('Dissolution potential [V]')


    # data = np.vstack((U_add,U_kink,U_edge,U_100,U_111)).T

    # np.savetxt(f'{metal}_HEA_dissolution.csv',data,delimiter=',',header='4,6,7,8,9',comments='cn:')


for ax in axes[0]:
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

for ax in axes[1]:
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))

    # ax.set_xlabel('Coordination Number')
    ax.set_xlabel('CN')

axes[1,3].yaxis.set_major_locator(MultipleLocator(1))


axes[0,1].set_xticks([3,4,5,6,7,8,9])
axes[0,0].set_ylabel('$\Delta E$ [eV]')
axes[1,0].set_ylabel('U [V]')


handles = []
labels = []
for metal in metals:
    handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12)) 
    labels.append('')
    handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12)) 
    labels.append('')
    handles.append(Line2D([0], [0], marker='o', color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=12)) 
    labels.append(metal)

handles.append(Line2D([0], [0], marker='o', color="w", label='',markerfacecolor='w', markersize=12))
labels.append('')

plt.tight_layout()
plt.subplots_adjust(top=0.9)

pos1 = axes[0,0].get_position()
pos2 = axes[0,4].get_position()

# mean_pos1 = (pos1.x0 + pos1.x1)/2
# mean_pos2 = (pos2.x0 + pos2.x1)/2

fig.legend(handles=handles[1:], labels=labels[1:],
           loc='outside upper center', ncol=len(handles), mode='expand',fontsize=16,bbox_to_anchor=(pos1.x0, .5, pos2.x1-pos1.x0, 0.5),fancybox=False)


plt.savefig('DFT_violin.png',dpi=600,bbox_inches='tight')
