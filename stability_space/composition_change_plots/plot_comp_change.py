import numpy as np
import matplotlib.pyplot as plt
from utilities.colors import alloy_color
from utilities.compositionspace_functions import prepare_triangle_plot, get_molar_fractions, molar_fractions_to_cartesians, prepare_tetrahedron_plot
from utilities import metals


mfs = get_molar_fractions(0.01,3)
xgrid = molar_fractions_to_cartesians(mfs)

data = np.loadtxt('grid_search/grid_data.csv',delimiter=',',skiprows=1)
bulk_comp_all = data[:,:8]
isc_all = data[:,8:16]
fsc_all = data[:,16:24]
diss_all = data[:,24:-1]
Sd_all = data[:,-1]

sum_ = np.sum(fsc_all,axis=1)
mask = (sum_<0.99) * (sum_>0)

fsc_all[mask] /=np.sum(fsc_all[mask],axis=1).reshape(-1,1)
diss_all[mask] /=np.sum(diss_all[mask],axis=1).reshape(-1,1)




quaternary_metals = ['Au','Pd','Pt','Cu']
system = ''.join(quaternary_metals)

metal_mask = [metal in quaternary_metals for metal in metals]

bulk_comp = bulk_comp_all[:,metal_mask]
isc = isc_all[:,metal_mask]
fsc = fsc_all[:,metal_mask]
diss = diss_all[:,metal_mask]

# Quarternary mask including only 12.5 at.% grid
q_mask = np.isclose(np.sum(bulk_comp,axis=1),1)*np.all(np.isclose(bulk_comp%0.125,0.0),axis=1)

bulk_comp = bulk_comp[q_mask]
isc = isc[q_mask]
fsc = fsc[q_mask]
diss = diss[q_mask]
Sd = Sd_all[q_mask]

# Change order of metals
sort_ids = [0,2,3,1]
bulk_comp = bulk_comp[:,sort_ids]
isc = isc[:,sort_ids]
fsc = fsc[:,sort_ids]
diss = diss[:,sort_ids]

bulk_colors = np.array([alloy_color(quaternary_metals,mf) for mf in bulk_comp])

# Stabiliy masks
stab_mask = Sd>0.0
unstab_mask = np.invert(stab_mask)


# Plot tetrahedron
ax = plt.figure(figsize=(8,8)).add_subplot(projection='3d')
ax = prepare_tetrahedron_plot(ax,quaternary_metals)

X = molar_fractions_to_cartesians(fsc[stab_mask]).T
X_unstab = molar_fractions_to_cartesians(bulk_comp[unstab_mask]).T

# Plot the final surface compositon with marker scaling with Sd
for i,x in enumerate(X.T):
    ax.scatter(x[0],x[1],x[2],color=bulk_colors[stab_mask][i],marker='o',s=50*Sd[stab_mask][i]+10,alpha=0.8,edgecolors='none')
# Plot Sd=0 as infulled markers
ax.scatter(X_unstab[0],X_unstab[1],X_unstab[2],facecolors='none',edgecolors=bulk_colors[unstab_mask],marker='o')

# Plot line connecting initial to final composition
X_isc = molar_fractions_to_cartesians(isc[stab_mask])
for i in range(X.shape[1]):
    ax.plot([X_isc[i,0],X[0,i]],[X_isc[i,1],X[1,i]],[X_isc[i,2],X[2,i]],c='k',alpha=0.08,zorder=0)

ax.view_init(12,-80)

plt.savefig(f'composition_change_plots/{system}_comp_change.svg',bbox_inches='tight')


# Ternary plot excluding Cu
sub_mask = np.isclose(fsc[:,3],0.0)
sub_mask_b = np.isclose(bulk_comp[:,3],0.0)

fig,ax = plt.subplots(figsize=(4,4))
ax = prepare_triangle_plot(ax,quaternary_metals[:3])

X = molar_fractions_to_cartesians(fsc[stab_mask*sub_mask,:3]).T
X_unstab = molar_fractions_to_cartesians(bulk_comp[unstab_mask*sub_mask_b,:3]).T

for i,x in enumerate(X.T):
    ax.scatter(x[0],x[1],color=bulk_colors[stab_mask*sub_mask][i],marker='o',s=60*Sd[stab_mask*sub_mask][i]+5,alpha=0.8,edgecolors='none')
ax.scatter(X_unstab[0],X_unstab[1],facecolors='none',edgecolors=bulk_colors[unstab_mask*sub_mask_b],marker='o')


plt.savefig(f'composition_change_plots/{system}_ternarysubspace_comp_change.png',dpi=600,bbox_inches='tight')    
