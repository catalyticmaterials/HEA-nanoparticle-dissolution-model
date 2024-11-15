from utilities.particle_dissolver import Dissolver
import numpy as np
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write
from ase.data import covalent_radii, atomic_numbers

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)

Sd = []
dummy = dissolver.dummy_particle(1925)
N_initial_111,_ = dissolver.get_composition_cn(dummy,9)

metals=['Au','Pd']





def get_min_Udiss(atoms):
    CN,nl = dissolver.get_coordination_numbers(atoms,True)
    symbols = np.array(atoms.get_chemical_symbols())
    ids = np.arange(len(atoms))
    surface_mask = CN<=9
    surface_ids = ids[surface_mask]
    surface_symbols = np.array(symbols)[surface_mask]
    features = dissolver.linear_neighborid_features(nl,surface_ids,symbols,CN)
    dE = dissolver.get_dE(surface_symbols,features)

    dissolution_potentials = dissolver.get_Udiss(surface_symbols,dE)

    return np.min(dissolution_potentials)



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

    write(f'Au_protection/pov/{name}.pov', atoms, colors=colors,rotation=rotation,radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    



# Clean Pd
Pd_particle = dissolver.make_particle([0.0,0.0,0.0,0.0,1,0.0,0.0,0.0],1925,return_particle=True)
povray_figure(Pd_particle,f'Pd','12x,-12y,-3z')
CNs=dissolver.get_coordination_numbers(Pd_particle)
min_Udiss = [get_min_Udiss(Pd_particle)]

# Galvanic replacement
for cn in (6,7,8):
    PdAu_particle = Pd_particle.copy()

    symbols = np.array(['Pd']*len(Pd_particle))
    symbols[CNs<=cn] = 'Au'
    PdAu_particle.set_chemical_symbols(symbols)
    

    povray_figure(PdAu_particle,f'PdAu_{cn}','12x,-12y,-3z')

    min_Udiss.append(get_min_Udiss(PdAu_particle))




print(min_Udiss)
