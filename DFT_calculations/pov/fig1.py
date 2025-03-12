import numpy as np
from ase.visualize import view
from utilities import metals
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write, Trajectory
from ase.data import covalent_radii, atomic_numbers
from utilities.particle_dissolver import Dissolver



colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}

def get_colors(atoms):
    colors = np.zeros((len(atoms),4))
    symbols = np.array(atoms.get_chemical_symbols())
    for metal in np.unique(symbols):
        mask = symbols==metal
        colors[mask,:3] = [colors_dict[metal]]*np.sum(mask)
    return colors

def povray_figure(atoms,name,rotation,uc,tp=None):
    # bbox_coord = np.max(np.abs(atoms.get_positions()))*1.4
    colors = get_colors(atoms)
    
    if tp is not None:
        colors[40,3] = tp

    # radii = r_Pt*len(atoms)
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
    filename = f'pov/{name}.pov'
    write(filename, atoms, colors=colors,rotation=rotation,bbox = [-3, -3, 15, 35], show_unit_cell=uc, povray_settings=povray_settings)
    
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

from ase.db import connect

db = connect('hea_dbs/T111_HEA.db')

atoms = db.get_atoms(40)

povray_figure(atoms,'T111_example_tp','-60x',uc=1,tp=0.75)
povray_figure(atoms,'T111_example','-60x',uc=1,tp=0.0)
povray_figure(atoms,'T111_example_diss','-60x',uc=1,tp=1.0)

del atoms[np.delete(np.arange(45),40)]

povray_figure(atoms,'Cu_target','-60x',uc=0)








