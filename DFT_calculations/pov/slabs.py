import numpy as np
from ase.visualize import view
from utilities import metals
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write, Trajectory
from ase.data import covalent_radii, atomic_numbers
from utilities.particle_dissolver import Dissolver
from ase.db import connect



colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}

def get_colors(atoms):
    colors = np.zeros((len(atoms),4))
    symbols = np.array(atoms.get_chemical_symbols())
    for metal in np.unique(symbols):
        mask = symbols==metal
        colors[mask,:3] = [colors_dict[metal]]*np.sum(mask)
    return colors

def povray_figure(atoms,name,rotation,bbox,colors):
    # bbox_coord = np.max(np.abs(atoms.get_positions()))*1.4
    # colors = get_colors(atoms)
    
    # if uc==1:
    #     colors[40,3] = 0.75

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
    filename = f'pov/slabs/{name}.pov'
    write(filename, atoms, colors=colors,rotation=rotation,bbox = bbox, show_unit_cell=0, povray_settings=povray_settings)
    
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




aoi_color = 0.2



# CN9
db = connect('hea_dbs/T111_HEA.db')

atoms = db.get_atoms(40)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[40,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN9','-70x',[-3,7,50,30],atom_colors_rep)


# CN3
db = connect('hea_dbs/T111_HEA.db')

atoms = db.get_atoms(38)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[45,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN3','-70x',[-3,7,50,32],atom_colors_rep)





# CN8
db = connect('hea_dbs/T100_HEA.db')

atoms = db.get_atoms(40)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[40,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN8','-70x,-30y,-3z',[-5,5,45,30],atom_colors_rep)


# CN4
db = connect('hea_dbs/T100_HEA.db')

atoms = db.get_atoms(38)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[45,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN4','-70x,-30y,-3z',[-5,5,45,30],atom_colors_rep)





# CN7
db = connect('hea_dbs/Edge_HEA.db')

atoms = db.get_atoms(39)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[1,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN7','-75x,-15y,-3z',[-5,5,45,30],atom_colors_rep)


# CN5
db = connect('hea_dbs/Edge_HEA.db')

atoms = db.get_atoms(40)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[45,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN5','-75x,-15y,-3z',[-5,5,45,30],atom_colors_rep)




# CN6
db = connect('hea_dbs/Kink_HEA.db')

atoms = db.get_atoms(19)

atom_colors = np.zeros((len(atoms),4))

atom_colors[:,:3] = np.array([colors_dict[symbol] for symbol in atoms.get_chemical_symbols()])
atom_colors[:,:3] += (1 - atom_colors[:,:3])*0.5

atom_colors[2,:3] = aoi_color

atoms_rep = atoms.repeat((4,4,1))
atom_colors_rep = np.vstack([atom_colors]*16)

povray_figure(atoms_rep,'CN6','63x,-22y,-162z',[-30,-10,25,25],atom_colors_rep)