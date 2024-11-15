from utilities.particle_dissolver import Dissolver
import numpy as np
from utilities import metals
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write
from ase.data import covalent_radii, atomic_numbers

import matplotlib.pyplot as plt

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
np.random.seed(2)
dissolver.make_particle(np.ones(8)/8,380)

atoms,diss,traj=dissolver.dissolve_atoms_batch(0.8,relax_func=dissolver.relax_particle_batch_cn,return_trajectory=True)



for i,poi in enumerate(traj):
    cn,nl = dissolver.get_coordination_numbers(poi,True)

    if np.any(cn==3):
        print(i)
        break




cn,nl = dissolver.get_coordination_numbers(poi,True)

surface_mask = (cn<=9)*(cn>=3)

symbols = np.array(poi.get_chemical_symbols())
surface_symbols = symbols[surface_mask]

features = dissolver.linear_neighborid_features(nl,np.where(surface_mask)[0],symbols,cn)


dE = dissolver.get_dE(surface_symbols,features)

Udiss = dissolver.get_Udiss(surface_symbols,dE) - 0.8







bbox_coord = np.max(np.abs(traj[0].get_positions()))*1.4


colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}



def get_metal_colors(atoms):
    colors = np.zeros((len(atoms),4))
    symbols = np.array(atoms.get_chemical_symbols())
    for metal in np.unique(symbols):
        mask = symbols==metal
        colors[mask,:3] = [colors_dict[metal]]*np.sum(mask)
    return colors


cn_colors_dict = {
                1: to_rgb('brown'),
                2: to_rgb('brown'),
                3: to_rgb('deeppink'),
                4: to_rgb('lightskyblue'),
                5: to_rgb('darkorchid'),
                6:to_rgb('darkorange'),
                7:to_rgb('blue'),
                8:to_rgb('yellow'),
                9:to_rgb('limegreen'),
                10:to_rgb('dimgrey'),
                11:to_rgb('dimgrey'),
                12:to_rgb('dimgrey')}

def get_cn_colors(atoms):
    colors = np.zeros((len(atoms),4))
    CN = dissolver.get_coordination_numbers(atoms)
    for cn in np.unique(CN):
        mask = CN==cn
        colors[mask,:3] = [cn_colors_dict[cn]]*np.sum(mask)
    return colors


r_Pt = [covalent_radii[atomic_numbers['Pt']]]

def povray_figure(atoms,colors,name):
    
    
    radii = r_Pt*len(atoms)
    povray_settings = {
        'display': False,  # Display while rendering
        'pause': True,  # Pause when done rendering (only if display)
        'transparent': True,  # Transparent background
        # 'canvas_width': 1500,  # Width of canvas in pixels
        'canvas_height': 1500,  # Height of canvas in pixels
        #'camera_dist': 50.,  # Distance from camera to front atom
        'image_plane': None,  # Distance from front atom to image plane
        #'camera_type': 'perspective',  # perspective, ultra_wide_angle
        #'point_lights': [[(0,0,0),'White'],[(-30,0,0),'White']],             # [[loc1, color1], [loc2, color2],...]
        #'point_lights': [[(0,20,0),'Gray50']],  
        'area_light': [(0, 0, 45),  # location
                    'White',       # color
                    1, 1, 3, 3,],  # width, height, Nlamps_x, Nlamps_y
        'background': 'Grey',        # color
        #'textures': None,  # Length of atoms list of texture names
        'transmittances': colors[:,3],
        'celllinewidth': 0.1,  # Radius of the cylinders representing the cell
        # 'image_width': 1000,
        # 'image_height': 1000 
       

    }

    filename = f'demonstrate_regressor/pov/{name}.pov'
    write(filename, atoms, colors=colors,rotation='12x,-12y,-3z',radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    
    # Costum light setting
    with open(filename, "r") as f:
        lines = f.readlines()
    lines[10:15] = ['light_source {<  -30.00,  30.00,   40.00> color Gray40 shadowless}\n',
                             'light_source {<  30.00,  30.00,   40.00> color Gray40 shadowless}\n',
                             'light_source {<  30.0,  -30.00,   40.00> color Gray40 shadowless}\n',
                             'light_source {<  -30.0,  -30.00,   40.00> color Gray40 shadowless}\n',
                             'light_source {<  0.0,  0.00,   40.00> color Gray25 shadowless}\n']
    with open(filename, "w") as f:
        for line in lines:
            f.write(line)
    


metal_colors = get_metal_colors(poi)
povray_figure(poi,metal_colors,'small_metalcolors')

import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
# make colorbar
# Create a ListedColormap from the colors
cmap = mcolors.ListedColormap([colors_dict[metal] for metal in metals])

# Define the boundaries for each color
bounds = np.linspace(0, 8, 9)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Plot the colorbar
fig, ax = plt.subplots(figsize=(0.15, 2))  # Adjust size as needed
cbar = plt.colorbar(
    ScalarMappable(cmap=cmap, norm=norm), cax=ax,
    ticks=np.arange(0.5, 8.5,1.0)  # Place ticks at the center of each color
)
cbar.ax.set_yticklabels(metals)  # Set custom labels
plt.minorticks_off()
plt.savefig('demonstrate_regressor/colorbars/metal_colorbar.svg',dpi=600,bbox_inches='tight')



cn_colors = get_cn_colors(poi)
povray_figure(poi,cn_colors,'small_cncolors')

# make colorbar
# Create a ListedColormap from the colors
cmap = mcolors.ListedColormap([cn_colors_dict[n] for n in np.arange(3,11)])

# Define the boundaries for each color
bounds = np.linspace(0, 8, 9)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Plot the colorbar
fig, ax = plt.subplots(figsize=(0.15, 2))  # Adjust size as needed
cbar = plt.colorbar(
    ScalarMappable(cmap=cmap, norm=norm), cax=ax,
    ticks=np.arange(0.5, 8.5,1.0)  # Place ticks at the center of each color
)
cn_labels = ['3','4','5','6','7','8','9','>9']
cbar.ax.set_yticklabels(cn_labels)  # Set custom labels
cbar.set_label(label='CN')
plt.minorticks_off()
plt.savefig('demonstrate_regressor/colorbars/cn_colorbar.svg',dpi=600,bbox_inches='tight')






max_dE = np.max(np.abs(dE))
max_Udiss = np.max(np.abs(Udiss))

import matplotlib.pyplot as plt
dE_norm = Normalize(-max_dE,max_dE)
cmap = plt.get_cmap('seismic_r')



dE_colors = np.zeros((len(poi),4))

dE_colors[cn>9,:3] = [to_rgb('dimgrey')[:3]]*np.sum(cn>9)
if np.any(cn<3):
    dE_colors[cn<3,:3] = [to_rgb('black')[:3]]*np.sum(cn<3)

dE_colors[surface_mask,:3] = cmap(dE_norm(dE))[:,:3] 

povray_figure(poi,dE_colors,'small_dEcolors')

# Coloarbar
fig,ax = plt.subplots(figsize=(0.15,2.4))
cbar=plt.colorbar(ScalarMappable(cmap=cmap,norm=dE_norm),cax=ax,orientation='vertical')
cbar.set_label(label=r'$\Delta E$ [eV]')
cbar.ax.minorticks_on()
plt.savefig('demonstrate_regressor/colorbars/dE_colorbar.svg',dpi=600,bbox_inches='tight')





Udiss_norm = Normalize(-max_Udiss,max_Udiss)

Udiss_colors = np.zeros((len(poi),4))

Udiss_colors[cn>9,:3] = [to_rgb('dimgrey')[:3]]*np.sum(cn>9)
if np.any(cn<3):
    Udiss_colors[cn<3,:3] = [to_rgb('black')[:3]]*np.sum(cn<3)

Udiss_colors[surface_mask,:3] = cmap(Udiss_norm(Udiss))[:,:3] 

povray_figure(poi,Udiss_colors,'small_Udisscolors')

# Coloarbar
fig,ax = plt.subplots(figsize=(0.15,2))
cbar=plt.colorbar(ScalarMappable(cmap=cmap,norm=Udiss_norm),cax=ax,orientation='vertical')
cbar.set_label(label=r'$U_{diss} - U$ [V]')
cbar.ax.minorticks_on()
plt.savefig('demonstrate_regressor/colorbars/Udiss_colorbar.svg',dpi=600,bbox_inches='tight')