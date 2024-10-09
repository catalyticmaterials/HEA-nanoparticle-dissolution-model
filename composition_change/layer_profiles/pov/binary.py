import numpy as np
from ase.visualize import view
from utilities import metals
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write, Trajectory
from ase.data import covalent_radii, atomic_numbers




colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}
colors_dict['Pd'] = (10/255,15/255,225/255)
colors_dict['Pt'] = (0/255,100/255,0/255)
colors_dict['Ir'] = (40/255,0,85/255)
colors_dict['Ru'] = (120/255,0,0/255)


def get_colors(atoms):
    colors = np.zeros((len(atoms),4))
    symbols = np.array(atoms.get_chemical_symbols())
    for metal in np.unique(symbols):
        mask = symbols==metal
        colors[mask,:3] = [colors_dict[metal]]*np.sum(mask)
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

    write(f'layer_profiles/pov/{name}.pov', atoms, colors=colors,rotation='0x,-30y,0z',radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    




# Make cross section
def cross_section(atoms_):
    atoms = atoms_.copy()
    pos = atoms.get_positions()
    mask = pos[:,0]>0.0
    del atoms[mask]
    return atoms



traj = Trajectory('layer_profiles/Pt75Cu25.traj')

initial = traj[0]
final = traj[-1]

initial_cs = cross_section(initial)
final_cs = cross_section(final)

povray_figure(initial_cs,'Pt75Cu25_i_cs')
povray_figure(final_cs,'Pt75Cu25_f_cs')






traj = Trajectory('layer_profiles/Pt25Cu75.traj')

initial = traj[0]
final = traj[-1]

initial_cs = cross_section(initial)
final_cs = cross_section(final)

povray_figure(initial_cs,'Pt25Cu75_i_cs')
povray_figure(final_cs,'Pt25Cu75_f_cs')






traj = Trajectory('layer_profiles/PtRu.traj')

initial = traj[0]
final = traj[-1]

initial_cs = cross_section(initial)
final_cs = cross_section(final)

povray_figure(initial_cs,'PtRu_i_cs')
povray_figure(final_cs,'PtRu_f_cs')




traj = Trajectory('layer_profiles/PdRu.traj')

initial = traj[0]
final = traj[-1]

initial_cs = cross_section(initial)
final_cs = cross_section(final)

povray_figure(initial_cs,'PdRu_i_cs')
povray_figure(final_cs,'PdRu_f_cs')