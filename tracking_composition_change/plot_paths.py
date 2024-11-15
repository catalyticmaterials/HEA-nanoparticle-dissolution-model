from utilities.compositionspace_functions import prepare_triangle_plot, molar_fractions_to_cartesians
import numpy as np
import matplotlib.pyplot as plt

n=10

data = np.loadtxt(f'initial_and_final_111_comps.csv',delimiter=',',skiprows=1)
traj_data = [np.loadtxt(f'trajectories/traj_111_comp_{i}.csv',delimiter=',',skiprows=1) for i in range(n)]

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


# Plot paths together
fig,ax = plt.subplots(figsize=(4,4))
ax = prepare_triangle_plot(ax,ternary_metals)
for i in range(n):
    traj_i = pseudo_mfs(traj_data[i][:,:8])
    diss_i = pseudo_mfs(traj_data[i][1:,8:])

    X_traj = molar_fractions_to_cartesians(traj_i)
    X_diss = molar_fractions_to_cartesians(diss_i)

    ax.plot(X_traj.T[0],X_traj.T[1],c='k',alpha=0.2,zorder=0)
    ax.scatter(X_traj.T[0,0],X_traj.T[1,0],c='tab:blue',edgecolor='blue',zorder=1)
    ax.scatter(X_traj.T[0,-1],X_traj.T[1,-1],c='tab:red',alpha=0.75,marker='x')


    ax.plot(X_diss.T[0],X_diss.T[1],c='k',alpha=0.2,zorder=0)
    ax.scatter(X_diss.T[0,0],X_diss.T[1,0],c='tab:green',edgecolor='green',zorder=1,marker='^')
    ax.scatter(X_diss.T[0,-1],X_diss.T[1,-1],c='tab:red',alpha=0.75,marker='x')


plt.savefig('paths.png',dpi=600,bbox_inches='tight')


# Plot paths individually
fig,axes = plt.subplots(ncols=2,nrows=5,figsize=(9,18))
axes = axes.flatten()

for i,ax in enumerate(axes):
    traj_i = pseudo_mfs(traj_data[i][:,:8])
    diss_i = pseudo_mfs(traj_data[i][1:,8:])

    X_traj = molar_fractions_to_cartesians(traj_i)
    X_diss = molar_fractions_to_cartesians(diss_i)

    ax = prepare_triangle_plot(ax,ternary_metals)

    ax.plot(X_traj.T[0],X_traj.T[1],c='k',alpha=0.2,zorder=0)
    ax.scatter(X_traj.T[0,0],X_traj.T[1,0],c='tab:blue',edgecolor='blue',zorder=1)
    ax.scatter(X_traj.T[0,-1],X_traj.T[1,-1],c='tab:red',alpha=0.75,marker='x')


    ax.plot(X_diss.T[0],X_diss.T[1],c='k',alpha=0.2,zorder=0)
    ax.scatter(X_diss.T[0,0],X_diss.T[1,0],c='tab:green',edgecolor='green',zorder=1,marker='^')
    ax.scatter(X_diss.T[0,-1],X_diss.T[1,-1],c='tab:red',alpha=0.75,marker='x')


plt.subplots_adjust(wspace=0.4,hspace=0.3)

plt.savefig('paths_individual.png',dpi=600,bbox_inches='tight')

