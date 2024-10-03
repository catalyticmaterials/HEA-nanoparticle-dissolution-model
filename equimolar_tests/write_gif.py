from utilities import metals
from utilities import metal_colors
from utilities.particle_dissolver import Dissolver
from ase.io import write, Trajectory
import numpy as np
from matplotlib.colors import to_rgb
from ase.data import covalent_radii, atomic_numbers
import imageio


# dissolver = Dissolver(metals,regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
# dissolver.make_particle([1/8]*8)

# Dissolver.dissolve_atoms(0.8,write_traj='equimolar_norelax.traj',relax_cn=False)
# dissolver.dissolve_atoms(0.8,traj_file='equimolar_relax.traj',relax_func=dissolver.relax_particle_batch_cn)


traj_relax = Trajectory('equimolar_relax.traj')
# traj_norelax = Trajectory('equimolar_norelax.traj')

# write('equimolar_relax.gif',traj_relax,interval=200)
# write('equimolar_norelax.gif',traj_norelax,interval=200)






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
        'transparent': False,  # Transparent background
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
        # 'transmittances': colors[:,3],
        'celllinewidth': 0.1,  # Radius of the cylinders representing the cell
        # 'image_width': 1000,
        # 'image_height': 1000
    }

    write(f'images/{name}.pov', atoms, colors=colors,rotation='12x,-12y,-3z',radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    





# for i,atoms in enumerate(traj_relax):

#     povray_figure(atoms,f'eqm_relax_{i}')


frames = [imageio.imread(f'images/eqm_relax_{i}.png') for i in range(len(traj_relax))]

imageio.mimsave('eqm_relax_gif.gif',frames,fps=2)