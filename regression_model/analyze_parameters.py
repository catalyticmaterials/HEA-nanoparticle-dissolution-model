import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from utilities import metals, U_standard, load_regressor
from utilities.colors import metal_colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ase.db import connect
from matplotlib.lines import Line2D


regressor = load_regressor('../utilities/AgAuCuIrPdPtRhRu_multilinear.regressor')


metal_data = np.loadtxt('../DFT_calculations/metals_data.csv',delimiter=',',skiprows=1,usecols=np.arange(16)).T

markers = ['P','p','H','s','o','^','X','d']

fig,ax = plt.subplots(figsize=(5,3.2))
cns = np.arange(3,10)
diss_cont = []
diffs=[]
E_diff=[]
all_parameters = []
for i,metal in enumerate(metals):
    # print(metal)
    params = regressor.parameters[metal]
    all_parameters.append(params)

    cn_params = params[:7]
    neighbor_params = params[7:]

    diss_cont.append(neighbor_params)
    
    cn_params_U = np.min(cn_params/U_standard[metal]['n'] + U_standard[metal]['U'],axis=0)

    metal_U = metal_data[i+8]

    ax.scatter(cns,cn_params_U, marker=markers[i],color=metal_colors[metal],label=metal)
    # ax.scatter(cns,cn_params, marker=markers[i],color=metal_colors[metal],label=metal)
    ax.scatter(cns,metal_U, marker=markers[i],alpha=0.6,color=metal_colors[metal])

    diffs.append(cn_params_U - metal_U)

    E_diff.append(cn_params - metal_data[i])


np.savetxt('parameters/parameters.csv',np.array(all_parameters).T,fmt='%1.2f',delimiter=',',header=','.join(metals))



ME_metal = np.mean(diffs,axis=1)
MAE_metal = np.mean(np.abs(diffs),axis=1)
for i, metal in enumerate(metals):
    print(metal,ME_metal[i],MAE_metal[i])


diffs = np.array(diffs).flatten()
ME = np.mean(diffs)
MAE = np.mean(np.abs(diffs))

ax.set_xlabel('Coordination Number')
ax.set_ylabel('Dissolution potential [V]')

ax.text(0.02,0.98,f'Parameter vs. metal DFT:\nME: {ME:1.3f} [V], MAE: {MAE:1.3f} [V]',va='top',transform=ax.transAxes)
ax_ins = ax.inset_axes((0.58,0.2,0.4,0.3))
ax_ins.hist(diffs,bins=8,color='k',alpha=0.6)
ax_ins.set_facecolor('none')
ax_ins.spines['right'].set_visible(False)
ax_ins.spines['top'].set_visible(False)
ax_ins.spines['left'].set_visible(False)
ax_ins.yaxis.tick_left()
ax_ins.tick_params(axis='y', left=False, labelleft=False)
ax_ins.tick_params(axis='x', labelsize='small')
ax_ins.set_xlabel('Parameters - $U_{metals}^{DFT}$',size='small')


h, l = ax.get_legend_handles_labels()

plt.tight_layout()
plt.subplots_adjust(top=0.875)
axpos = ax.get_position()

fig.legend(handles=h, labels=l,
           loc='outside upper center', ncols=8,mode='expand',bbox_to_anchor=(0.02, .5, 0.96, 0.5),fancybox=False)


plt.savefig('parameters/parameters_vs_dft.png',dpi=600,bbox_inches='tight')
plt.close()





parameters = np.array([regressor.parameters[metal] for metal in metals])

fig,ax = plt.subplots(figsize=(6,3))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.025)

# Determine the maximum absolute value in the data for symmetric color scaling
max_abs_value = np.max(np.abs(parameters))

# Create a heatmap with seismic colormap centered around 0
heatmap = ax.imshow(parameters,interpolation='none',cmap='seismic_r',vmin=-max_abs_value, vmax=max_abs_value)

cbar = plt.colorbar(heatmap,cax)
cbar.set_label('Parameters [eV]', rotation=270, labelpad=15)

ax.text(3,9,'CN',ha='center')
ax.text(10.5,9,'Neighbor',ha='center')
ax.set_ylabel('Target atom')

ax.set_xticks(np.arange(15),[3,4,5,6,7,8,9]+metals)
ax.set_yticks(np.arange(len(metals)),labels=metals)

ax.axvline(6.5,c='k')

plt.savefig('parameters/parameters_heatmap.png',dpi=600,bbox_inches='tight')
plt.close()







with connect('../DFT_calculations/metal_dbs/bulks.db') as db:

    E_bulk = np.array([row.energy for row in db.select()])


surfaces = ['T111','T100','edge','kink']
E_surfaces = {}
with connect('../DFT_calculations/metals_dbs/metals_dissolution_out.db') as db:
    for surface in surfaces:
        E_surface = []
        for i,row in enumerate(db.select(surface=surface,defect='none')):
            E_total = row.energy

            atoms = db.get_atoms(row.id)
            n=len(atoms)

            a,b = atoms.cell.cellpar()[:2]
            A = a*b

            E_surface.append((E_total - n*E_bulk[i])/(2*A))
        E_surfaces[surface] = np.array(E_surface)




fig,ax = plt.subplots(figsize=(5,3))
E_surface_rels = []

for i,metal in enumerate(metals):

    # E_bulk_rel = E_bulk - E_bulk[i]
    E_surface_rel = E_surfaces['T111'] - E_surfaces['T111'][i]
    E_surface_rels.append(E_surface_rel)

    for j, marker in enumerate(markers):
        ax.scatter(E_surface_rel[j],diss_cont[i,j],color=metal_colors[metal],marker=markers[j])

    
h=[]
l=[]
for metal,marker in zip(metals,markers):
    h.append(Line2D([0], [0], marker=marker, color="w", label=metal,markerfacecolor=metal_colors[metal], markersize=8))
    l.append(metal)




ax.set_xlabel(r'$\gamma_{neighbor}^{(111)} - \gamma_{target}^{(111)}$ [eV]')
ax.set_ylabel('Energy Perturbation [eV]')



ax.xaxis.set_major_locator(MultipleLocator(0.05))
ax.xaxis.set_minor_locator(MultipleLocator(0.01))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))


plt.tight_layout()
plt.subplots_adjust(top=0.9)
axpos = ax.get_position()
# plt.legend(handles, labels, ncol=6,loc='upper left')
fig.legend(handles=h, labels=l,
           loc='outside upper center', ncols=8,fontsize=8,markerscale=1,mode='expand',bbox_to_anchor=(0.02, .5, 0.96, 0.5),fancybox=False)


plt.savefig('parameters/diss_cont_corr.png',dpi=600,bbox_inches='tight')