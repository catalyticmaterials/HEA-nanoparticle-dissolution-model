from utilities.particle_dissolver import Dissolver
from ase.data import covalent_radii, atomic_numbers
from ase.io import write
import numpy as np
from matplotlib.colors import to_rgb


colors_dict = {6:to_rgb('darkorange'),
               7:to_rgb('blue'),
               8:to_rgb('yellow'),
               9:to_rgb('limegreen'),
               12:to_rgb('dimgrey')}



dissolver = Dissolver()
particle = dissolver.dummy_particle()

def get_colors(atoms):
    colors = np.zeros((len(atoms),4))
    CN = dissolver.get_coordination_numbers(atoms)
    for cn in np.unique(CN):
        mask = CN==cn
        colors[mask,:3] = [colors_dict[cn]]*np.sum(mask)
    return colors

r_Pt = [covalent_radii[atomic_numbers['Pt']]]

def povray_figure(atoms,name):
    bbox_coord = np.max(np.abs(atoms.get_positions()))*1.4
    colors = get_colors(atoms)
    
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
        'point_lights': [],             # [[loc1, color1], [loc2, color2],...]
        'area_light': [(2., 3., 40.),  # location
                    'White',       # color
                    .7, .7, 3, 3],  # width, height, Nlamps_x, Nlamps_y
        'background': 'White',        # color
        #'textures': None,  # Length of atoms list of texture names
        'transmittances': colors[:,3],
        'celllinewidth': 0.1,  # Radius of the cylinders representing the cell
        # 'image_width': 1000,
        # 'image_height': 1000
    }
    # rotation='12x,-12y,-3z'
    write(f'{name}.pov', atoms, colors=colors,radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    

povray_figure(particle,'CN_particle')