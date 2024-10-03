from utilities.particle_dissolver import Dissolver
import numpy as np
from tqdm import tqdm
from ase.visualize import view
import matplotlib.pyplot as plt
from utilities.colors import metal_colors
from matplotlib.colors import to_rgb
from ase.io import write
from ase.data import covalent_radii, atomic_numbers

dissolver = Dissolver(regressor='AgAuCuIrPdPtRhRu_multilinear',c_metals=1e-6)
n=2
Sd = []
dummy = dissolver.dummy_particle(700)
N_initial_111,_ = dissolver.get_composition_cn(dummy,9)
# view(dummy)
# stop
# f_swap = 0.05
# n_swap = int(len(dummy)*f_swap)
# print(n_swap)
# view(dummy)
metals=['Au','Pd']





def get_min_Udiss(atoms):
    CN,nl = dissolver.get_coordination_numbers(atoms,True)
    symbols = np.array(atoms.get_chemical_symbols())
    ids = np.arange(len(atoms))
    surface_mask = CN<=9
    surface_ids = ids[surface_mask]
    surface_symbols = np.array(symbols)[surface_mask]
    features = dissolver.linear_neighborid_features(nl,surface_ids,symbols,CN)
    # print(np.unique(features,axis=0))
    # dissolution_potentials = dissolver.get_dissolution_potential(surface_symbols,features)
    dE = dissolver.get_dE(surface_symbols,features)
    # print(np.unique(dE,return_counts=True))
    dissolution_potentials = dissolver.get_Udiss(surface_symbols,dE)
    # print(np.unique(dissolution_potentials,return_counts=True))
    # stop
    return np.min(dissolution_potentials)


Pd_particle = dissolver.make_particle([0.0,0.0,0.0,0.0,1,0.0,0.0,0.0],1300,return_particle=True)
CNs=dissolver.get_coordination_numbers(Pd_particle)
PdAu_particle = Pd_particle.copy()
symbols = np.array(['Pd']*len(Pd_particle))
symbols[CNs<=7] = 'Au'
PdAu_particle.set_chemical_symbols(symbols)
dissolver.particle = PdAu_particle.copy()
# atoms,_=dissolver.dissolve_atoms_mc(1.0,relax_cn=True,traj_file='PdAu_7_1V.traj')
print(get_min_Udiss(PdAu_particle))
atoms,_=dissolver.dissolve_atoms(0.95,relax_func=dissolver.relax_particle_batch_cn,traj_file='PdAu_7_095V_batch.traj')
atoms,_=dissolver.dissolve_atoms_it(0.95,relax_func=dissolver.relax_particle_single_cn,traj_file='PdAu_7_095V_it.traj')
atoms,_=dissolver.dissolve_atoms(1.0,relax_func=dissolver.relax_particle_batch_cn,traj_file='PdAu_7_1V_batch.traj')
atoms,_=dissolver.dissolve_atoms_it(1.0,relax_func=dissolver.relax_particle_single_cn,traj_file='PdAu_7_1V_it.traj')




stop


def dissolve(particle,U):
    dissolver.particle = particle
    atoms,_=dissolver.dissolve_atoms(U,relax_cn=True)

    N_final_111, final_111_comp = dissolver.get_composition_cn(atoms,9)
    # view(atoms)
    return N_final_111/N_initial_111, np.array(final_111_comp)[[1,4]], atoms



colors_dict = {metal:to_rgb(metal_colors[metal]) for metal in metals}
# colors_dict['Pd'] = (10/255,15/255,225/255)


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

    write(f'pov/{name}.pov', atoms, colors=colors,rotation=rotation,radii=radii,
      bbox = [-bbox_coord, -bbox_coord, bbox_coord, bbox_coord], show_unit_cell=0, povray_settings=povray_settings)
    



# Clean Pd
Pd_particle = dissolver.make_particle([0.0,0.0,0.0,0.0,1,0.0,0.0,0.0],1300,return_particle=True)
# povray_figure(Pd_particle,f'Pd','12x,-12y,-3z')
CNs=dissolver.get_coordination_numbers(Pd_particle)
min_Udiss = [get_min_Udiss(Pd_particle)]
Sds = [dissolve(Pd_particle,min_Udiss[-1]+0.1)[0]]
for cn in (6,7,8):
    PdAu_particle = Pd_particle.copy()

    symbols = np.array(['Pd']*len(Pd_particle))
    symbols[CNs<=cn] = 'Au'
    PdAu_particle.set_chemical_symbols(symbols)
    

    # povray_figure(PdAu_particle,f'PdAu_{cn}','12x,-12y,-3z')

    min_Udiss.append(get_min_Udiss(PdAu_particle))
    # view(PdAu_particle)
    Sd, fsc, diss_particle = dissolve(PdAu_particle,min_Udiss[-1]+0.1)

    Sds.append(Sd)
    
    # povray_figure(diss_particle,f'PdAu_{cn}_diss','12x,-12y,-3z')


    


# fig,ax = plt.subplots(figsize=(4,3))
# ax2 = ax.twinx()
# ax.scatter(range(4),min_Udiss)
# ax2.scatter(range(4),Sds,c='tab:orange')
print(Sds)
print(min_Udiss)
# plt.show()

# initial_particle = dissolver.make_particle([0.0,0.0,0.0,0.0,1,0.0,0.0,0.0],1300,return_particle=True)

# CN,nl = dissolver.get_coordination_numbers(initial_particle,True)
# symbols = np.array(initial_particle.get_chemical_symbols())
# ids = np.arange(len(initial_particle))
# surface_mask = CN<=9
# surface_ids = ids[surface_mask]
# surface_symbols = np.array(symbols)[surface_mask]
# features = dissolver.linear_neighborid_features(nl,surface_ids,symbols,CN)
# dissolution_potentials = dissolver.get_dissolution_potential(surface_symbols,features)
# print(np.unique(dissolution_potentials))
# # swap_ids = surface_ids[np.argsort(dissolution_potentials)[:n_swap]]
# swap_ids = CN==6
# n_swap = np.sum(swap_ids)
# f_swap = n_swap/len(initial_particle)
# symbols[swap_ids] = 'Au'
# initial_particle.set_chemical_symbols(symbols)
# view(dissolver.particle)

# atoms,_=dissolver.dissolve_atoms(1.0,relax_cn=True)

# N_final_111, final_111_comp = dissolver.get_composition_cn(atoms,9)

# Sd = N_final_111/N_initial_111


# print(Sd)
# print(final_111_comp)



# CN,nl = dissolver.get_coordination_numbers(atoms,True)
# symbols = np.array(atoms.get_chemical_symbols())
# ids = np.arange(len(atoms))
# surface_mask = CN<=9
# surface_ids = ids[surface_mask]
# surface_symbols = np.array(symbols)[surface_mask]
# features = dissolver.linear_neighborid_features(nl,surface_ids,symbols,CN)
# dissolution_potentials = dissolver.get_dissolution_potential(surface_symbols,features)
# print(np.unique(dissolution_potentials))
# view(atoms)