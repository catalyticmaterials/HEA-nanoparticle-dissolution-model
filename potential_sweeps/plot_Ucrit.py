from utilities.compositionspace_functions import make_ternary_plot, molar_fractions_to_cartesians, get_molar_fractions, truncate_colormap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from utilities import U_diss


data = []


for metals in (['Pd','Pt','Ru'],['Au','Cu','Pt'],['Au','Cu','Pd']):

    system = ''.join(metals)
    Ucs = np.loadtxt(f'{system}_potential_sweep.csv',delimiter=',',skiprows=1)[:,3:]
    U=np.max(Ucs,axis=1)

    data.append(U)


grid = molar_fractions_to_cartesians(get_molar_fractions(0.05,3)).T

Ucmax = np.max(data)
Ucmin = np.min(data)


for Uc,metals in zip(data,(['Pd','Pt','Ru'],['Au','Cu','Pt'],['Au','Cu','Pd'])):
    system = ''.join(metals)
    
    make_ternary_plot(grid,Uc,metals,vmin=Ucmin,vmax=Ucmax,colorbar=False,contour_levels=30,colormap='magma',maxval=1.0,minval=0.0)

    plt.savefig(f'{system}_Uc.png',dpi=600,bbox_inches='tight')

from matplotlib.colors import Normalize
cmap=truncate_colormap(plt.get_cmap('magma'),0.0,1.0)
fig,ax = plt.subplots(figsize=(0.2,3))
cbar=plt.colorbar(ScalarMappable(cmap=cmap,norm=Normalize(Ucmin,Ucmax)),cax=ax,orientation='vertical')
cbar.set_label(label='$U_c$',rotation='horizontal')
plt.savefig('Uc_colorbar.png',dpi=600,bbox_inches='tight')