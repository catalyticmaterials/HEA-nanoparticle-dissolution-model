import numpy as np
from utilities import metals
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write, Trajectory
from ase.data import covalent_radii, atomic_numbers
from utilities.particle_dissolver import Dissolver


particle = Dissolver.dummy_particle(1925)
bbox_coord = np.max(np.abs(particle.get_positions()))*1.4


colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}



def get_colors(atoms):
    colors = np.zeros((len(atoms),4))
    symbols = np.array(atoms.get_chemical_symbols())
    for metal in np.unique(symbols):
        mask = symbols==metal
        colors[mask,:3] = [colors_dict[metal]]*np.sum(mask)
    return colors

r_Pt = [covalent_radii[atomic_numbers['Pt']]]

def povray_figure(atoms,name,rotation):
    # bbox_coord = np.max(np.abs(atoms.get_positions()))*1.4
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

    filename = f'hea/pov/{name}.pov'
    write(filename, atoms, colors=colors,rotation=rotation,radii=radii,
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
    




# Make cross section
def cross_section(atoms_):
    atoms = atoms_.copy()
    pos = atoms.get_positions()
    mask = pos[:,0]>0.0
    del atoms[mask]
    return atoms





for i in range(10):
    traj = Trajectory(f'../tracking_composition_change/trajectories/traj_{i}.traj')


    initial = traj[0]
    final = traj[-1]

    initial_cs = cross_section(initial)
    final_cs = cross_section(final)

    povray_figure(initial_cs,f'eqm_i{i}_cs','0x,-90y,0z')
    povray_figure(final_cs,f'eqm_f{i}_cs','0x,-90y,0z')

    povray_figure(initial,f'eqm_i{i}','12x,-12y,-3z')
    povray_figure(final,f'eqm_f{i}','12x,-12y,-3z')


    if i==2:
        povray_figure(traj[1],f'eqm_2_1','12x,-12y,-3z')
        povray_figure(traj[2],f'eqm_2_2','12x,-12y,-3z')
        povray_figure(traj[3],f'eqm_2_3','12x,-12y,-3z')
        povray_figure(traj[4],f'eqm_2_4','12x,-12y,-3z')