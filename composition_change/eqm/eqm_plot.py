from utilities.compositionspace_functions import prepare_triangle_plot, molar_fractions_to_cartesians, get_molar_fractions, make_ternary_plot, truncate_colormap
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal


data = np.loadtxt(f'eqm/initial_and_final_111_comps.csv',delimiter=',',skiprows=1)

Sd = data[:,-1]


mf_initial = data[:,:8]
mf_final = data[:,8:16]
mf_diss = data[:,16:24]

def pseudo_mfs(mfs):
    return np.hstack((np.sum(mfs[:,[1,3,5]],axis=1).reshape(-1,1),np.sum(mfs[:,[0,4]],axis=1).reshape(-1,1),np.sum(mfs[:,[2,6,7]],axis=1).reshape(-1,1)))

pseudo_initial_mf = pseudo_mfs(mf_initial)
pseudo_final_mf = pseudo_mfs(mf_final)
pseudo_diss_mf = pseudo_mfs(mf_diss)

ternary_metals = ['AuIrPt','AgPd','CuRhRu']

x_initials = molar_fractions_to_cartesians(pseudo_initial_mf)
x_finals = molar_fractions_to_cartesians(pseudo_final_mf)
x_diss = molar_fractions_to_cartesians(pseudo_diss_mf)

cmap = truncate_colormap(plt.get_cmap('viridis'),minval=0.2)
fig,ax = plt.subplots(figsize=(4,4))
hist,xedges,yedges,image=ax.hist2d(x_diss[:,0],x_diss[:,1],range=[[0.0,1.0],[0.0,1.0]],bins=100,cmap=cmap)
ax.fill(np.array([0.0,0.0,1/np.sqrt(3)])-0.005,[0.0,1.0,1.0],color='white')
ax.fill(np.array([1.0,1.0,1-1/np.sqrt(3)])+0.005,[0.0,1.0,1.0],color='white')
ax = prepare_triangle_plot(ax,ternary_metals,False)

plt.savefig(f'eqm/diss_comps.png',dpi=600,bbox_inches='tight')


fig,ax = plt.subplots(figsize=(4,4))
ax = prepare_triangle_plot(ax,ternary_metals)
for i in range(10):
    traj_i = pseudo_mfs(np.loadtxt(f'eqm/trajectories/traj_111_comp_{i}.csv',delimiter=',',skiprows=1))

    X_traj = molar_fractions_to_cartesians(traj_i)

    ax.plot(X_traj.T[0],X_traj.T[1],c='k',alpha=0.2,zorder=0)

ax.scatter(x_initials.T[0],x_initials.T[1],c='tab:blue',edgecolor='blue',zorder=1)
ax.scatter(x_finals.T[0],x_finals.T[1],c='tab:red',alpha=0.75,marker='x')    

plt.savefig(f'eqm/paths.png',dpi=600,bbox_inches='tight')


fig,axes = plt.subplots(ncols=2,nrows=5,figsize=(9,18))
axes = axes.flatten()

for i,ax in enumerate(axes):
    traj_i = pseudo_mfs(np.loadtxt(f'eqm/trajectories/traj_111_comp_{i}.csv',delimiter=',',skiprows=1))

    X_traj = molar_fractions_to_cartesians(traj_i)

    prepare_triangle_plot(ax,ternary_metals)

    ax.plot(X_traj.T[0],X_traj.T[1],c='k',alpha=0.2,zorder=0)

    ax.scatter(x_initials.T[0,i],x_initials.T[1,i],c='tab:blue',edgecolor='blue',zorder=1)
    ax.scatter(x_finals.T[0,i],x_finals.T[1,i],c='tab:red',alpha=0.75,marker='x')    


plt.subplots_adjust(wspace=0.4,hspace=0.3)

plt.savefig(f'eqm/paths_sep.png',dpi=600,bbox_inches='tight')


mfs = get_molar_fractions(0.01,3)
x_grid = molar_fractions_to_cartesians(mfs)



# mean_initial = np.mean(x_initials,axis=0)
# cov_initial = np.cov(x_initials.T)

# mv_initial = multivariate_normal(mean_initial,cov=cov_initial)
# z_initial = mv_initial.pdf(x_grid)


# fig,ax = make_ternary_plot(x_grid.T,z_initial,ternary_metals,colorbar='True',cbar_label='pdf initial')

# plt.savefig(f'eqm_ternary/pdf_initial.png',dpi=600,bbox_inches='tight')



# mean_final = np.mean(x_finals,axis=0)
# cov_final = np.cov(x_finals.T)

# mv_final = multivariate_normal(mean_final,cov=cov_final,allow_singular=True)
# z_final = mv_final.pdf(x_grid)


# fig,ax = make_ternary_plot(x_grid.T,z_final,ternary_metals,colorbar='True',cbar_label='pdf initial')

# plt.savefig(f'eqm_ternary/pdf_final.png',dpi=600,bbox_inches='tight')


cmap = truncate_colormap(plt.get_cmap('viridis'),minval=0.2)
fig,ax = plt.subplots(figsize=(4,4))
ax.tricontourf(x_grid[:,0],x_grid[:,1], np.zeros(len(x_grid)),cmap=cmap,vmax=1)
hist,xedges,yedges,image=ax.hist2d(x_initials[:,0],x_initials[:,1],range=[[0.27,0.73],[0,0.46]],bins=20,cmap=cmap)
# hist,xedges,yedges,image=ax.hist2d(initials[:,0],initials[:,1],density=True,range=[[0.25,0.75],[0,0.5]],bins=25,cmap=truncate_colormap(plt.get_cmap('hot_r'),maxval=0.75))
ax = prepare_triangle_plot(ax,ternary_metals)
plt.colorbar(image,ax=ax,shrink=0.5,anchor=(0.0,0.85),label='Frequency')
plt.savefig(f'eqm/2dhist_initial.png',dpi=600,bbox_inches='tight')


fig,ax = plt.subplots(figsize=(4,4))
ax.tricontourf(x_grid[:,0],x_grid[:,1], np.zeros(len(x_grid)),cmap=cmap,vmax=1)
hist,xedges,yedges,image=ax.hist2d(x_finals[:,0],x_finals[:,1],range=[[0.27,0.73],[0.,0.46]],bins=23,cmap=cmap)
ax = prepare_triangle_plot(ax,ternary_metals)
plt.colorbar(image,ax=ax,shrink=0.5,anchor=(0.0,0.85),label='Frequency')
plt.savefig(f'eqm/2dhist_final.png',dpi=600,bbox_inches='tight')