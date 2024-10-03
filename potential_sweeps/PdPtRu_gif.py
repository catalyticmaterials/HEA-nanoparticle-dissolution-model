from utilities.compositionspace_functions import molar_fractions_to_cartesians, make_ternary_plot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.widgets import Slider
from matplotlib.cm import ScalarMappable
import imageio
import numpy as np


potentials = np.arange(0.4,1.3,0.01)

data = np.loadtxt('PdPtRu_potential_sweep_fd.csv',delimiter=',',skiprows=1)
# U_crit = np.loadtxt('PdPtRu_potential_sweep.csv',delimiter=',',skiprows=1,usecols=(3))

grid = molar_fractions_to_cartesians(data[:,:3]).T

fd_data = data[:,3:]



U_crit = []
for i in range(len(data)):
    fds = fd_data[i]
    U_crit_idx = np.where(fds==0)[0][0]
    U_crit.append(potentials[U_crit_idx])

U_max = np.max(U_crit)

fig,(ax1,ax2) = plt.subplots(ncols=2,figsize=(8,4),dpi=600)
ax_slider = plt.axes([0.1,0.05,0.8,0.03])
# cmap = plt.get_cmap('viridis')

ax2 = make_ternary_plot(grid,U_crit,['Pd','Pt','Ru'],ax=ax2,vmin=np.min(U_crit),vmax=np.max(U_crit),colorbar=True,contour_levels=30,cbar_label='$U_c$',cbar_ticks=[np.min(U_crit),0.6,0.7,0.8,0.9,1.0,1.1,np.max(U_crit)],colormap='magma',minval=0.0,maxval=0.95)

plt.colorbar(ScalarMappable(cmap='coolwarm_r'),ax=ax1,shrink=0.5,anchor=(0.0,0.85),ticks=np.arange(0,1.1,0.2),label='$f_d$')

frames=[]
for i,U in enumerate(potentials):

    fds = fd_data[:,i]

    ax1 = make_ternary_plot(grid,fds,['Pd','Pt','Ru'],ax=ax1,vmax=1,vmin=0,colorbar=False,colormap='coolwarm_r',minval=0.0,maxval=1.0)
 
    Slider(ax_slider,'U [V]',np.min(potentials),U_max+0.01,valinit=U)


    FigureCanvas(fig).draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    frames.append(image)


    ax1.clear()
    ax_slider.clear()
    



# save gif
imageio.mimsave('PdPtRu_potential_sweep.gif', frames, fps=3)