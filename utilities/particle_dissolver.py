import numpy as np
from ase.cluster import wulff_construction
from ase.neighborlist import NeighborList
from copy import deepcopy
from ase.io import Trajectory
from iteround import saferound
from . import U_standard, metals, load_regressor
import pkg_resources
from collections import deque

kB=8.617333262e-5
T=298.15


class Dissolver():
    def __init__(self,metals=metals,U_dict=U_standard,c_metals=1, feature_generator = None, regressor=None) -> None:
        self.metals = metals
        self.particle = None
        
        if feature_generator is None:
            self.feature_generator = self.linear_neighborid_features
        else: 
            self.feature_generator = feature_generator
        
        if regressor is not None:
            # Load regressor if giveb as name of regressor in utilities package
            if type(regressor)==str:
                try:
                    self.regressor = load_regressor(pkg_resources.resource_filename('utilities',f'{regressor}.regressor'))
                except:
                    raise ValueError('regressor not found')
            else:
                self.regressor = regressor
        
        self.nb_dist = 3.0
        self.min_cn = 3
        self.max_cn = 9
        self.n_cnfeat = 7

        if c_metals!=1:
            # Adjust Standard reduction potentials by ion concentration
            for metal in metals:
                U_dict[metal]['U'] = U_dict[metal]['U'] + kB*T/U_dict[metal]['n'] * np.log(c_metals)

        self.U_dict = U_dict
        
        

    def dummy_particle(self,n_atoms=1289,surfaces = [(1, 1, 1),(1,0,0)],energies=[1.,1.]):
        return wulff_construction('Pt', surfaces, energies, n_atoms,'fcc',rounding='closest')

    def make_particle(self,composition,n_atoms=1289,surfaces = [(1, 1, 1),(1,0,0)],energies=[1.,1.], return_particle=False):
        # Make particle
        particle = wulff_construction('Pt', surfaces, energies, n_atoms,'fcc',rounding='closest')
        # Draw elements randomly according to composition
        symbols = np.random.choice(self.metals,size=len(particle),p=composition)
        # Set symbols
        particle.set_chemical_symbols(symbols)

        self.particle = particle
        if return_particle:
            return particle
    
    def set_regressor(self,regressor):
        self.regressor = regressor
    
    
    def get_coordination_numbers(self,atoms, return_neighborlist=False):
        # make ASE neighbor list
        cutoff = [self.nb_dist/2]*len(atoms)
        nl = NeighborList(cutoff, self_interaction=False,bothways=True)
        nl.update(atoms)

        # get coordination numbers from len of neighbor IDs for each atom
        coordination_numbers = np.array([len(nl.get_neighbors(i)[0]) for i in range(len(atoms))],dtype=int)
        
        if return_neighborlist:
            return coordination_numbers, nl
        
        return coordination_numbers

    def linear_neighborid_features(self,nl,ids,symbols,coordination_numbers):

        neighbor_features = []
        for i in ids:
            # Get neighbor IDs and symbols
            neighbor_ids = nl.get_neighbors(i)[0]
            neighbor_symbols = symbols[neighbor_ids]
            # Get neighbor features as number of each metal
            neighbor_features.append([np.sum(neighbor_symbols==metal) for metal in self.metals])

        # Get index in CN feature by subtracting min CN from CNs
        coordination_ids = coordination_numbers[ids]-self.min_cn
        # Make CN features from Identity matrix
        coordination_features = np.eye(self.n_cnfeat,dtype=int)[coordination_ids]
        assert np.all(coordination_numbers[ids]>=3) and np.all(coordination_numbers[ids]<=9), f'{coordination_numbers[ids]},{ids}'
        # combine CN and neighbor features
        features=np.hstack((coordination_features, np.array(neighbor_features)/coordination_numbers[ids].reshape(-1,1)))

        return features
    
            

    def get_dissolution_potential(self,symbols,features):
        # initiate array
        U_diss = np.zeros(len(features))
        # Calculate U_diss by each target metal
        for metal in self.metals:
            mask = symbols==metal
            # U_diss = dE/ne + U_metal
            U_diss[mask] = np.min(self.regressor.predict(metal,features[mask])/self.U_dict[metal]['n'] + self.U_dict[metal]['U'],axis=0)
        return U_diss
                

    def get_composition(self,symbols,precision=2):
        # Count the number of each metal
        composition = [list(symbols).count(metal) for metal in self.metals]
        # Number of atoms
        N = len(symbols)
        if N==0:
            return N,np.zeros(len(self.metals))
        # Get composition by normalizing by N. Use saferound to asure summation to 1
        composition = saferound(np.array(composition)/N,precision)
        return N, composition

    def get_composition_cn(self,atoms,cn,precision=2):
        # get CN of all atoms
        cns = self.get_coordination_numbers(atoms)
        # Get symbols of specified CN by applying mask
        mask_cn = cns==cn
        symbols = np.array(atoms.get_chemical_symbols())[mask_cn]
        # Get composition
        return self.get_composition(symbols,precision)
    
    def composition_by_layer(self,atoms_,min_atoms=1,n_layers=None,precision=3):
        atoms = atoms_.copy()
        layer_compositions = []
        n_atoms_layer = []
        while len(atoms)>min_atoms:
            
            CN = self.get_coordination_numbers(atoms)

            surface_mask = CN<10

            symbols = np.array(atoms.get_chemical_symbols())
            layer_compositions.append(self.get_composition(symbols[surface_mask],precision=precision)[1])
            n_atoms_layer.append(np.sum(surface_mask))

            del atoms[surface_mask]

            if n_layers is not None and len(n_atoms_layer)==(n_layers -1):
                symbols = np.array(atoms.get_chemical_symbols())
                layer_compositions.append(self.get_composition(symbols,precision=precision)[1])
                n_atoms_layer.append(len(atoms))
                break
        else:
            layer_compositions = [[0]*len(self.metals)]*n_layers
            n_atoms_layer = [0]*n_layers

        return np.array(layer_compositions), np.array(n_atoms_layer)
        



    def update_cn(self,atoms,ids,nl,cn):
        nl.update(atoms)
        cn[ids] = [len(nl.get_neighbors(id_)[0]) for id_ in ids]

        return cn,nl

    def update_dE(self,dEs,ids,cn,nl,symbols):

        bulk_mask = cn[ids]>9
        dEs[ids[bulk_mask]] = np.inf 
        under_coord_mask = cn[ids]<3
        dEs[ids[under_coord_mask]] = -np.inf 
        
        ids = ids[np.invert(bulk_mask+under_coord_mask)]

        new_features = self.feature_generator(nl,ids,symbols,cn)

        for metal in self.metals:
            mask = symbols[ids]==metal
            dEs[ids[mask]] = self.regressor.predict(metal,new_features[mask])

        return dEs
    

    def get_cn_improvements(self,atoms,free_positions,CNs,atom_ids):
        # Get best improvement for each unoccupied position
        relax_to_id = []
        relax_atom_id = []
        improvements = []
        neighbors = []
        # Loop through free positions
        for i,free_position in enumerate(free_positions):
            
            dists = np.linalg.norm(atoms.positions - free_position,axis=1)

            nb_mask = dists<self.nb_dist

            n_nb = np.sum(nb_mask)

            if n_nb==0:
                relax_atom_id.append(None)
                relax_to_id.append(i)
                improvements.append(-100)
                neighbors.append([])
                continue


            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            # Change in CN
            d_CN = (n_nb - 1) - CNs[neighbor_ids]
           
            # Get improvement from moving a neighbor into the free position
            best_improvement = np.max(d_CN)

            # Move a neighbor into free position if it improves CN
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_CN)==best_improvement])
    
            # Save the neighbor to move, to where and its improvement
            relax_atom_id.append(idx)
            relax_to_id.append(i)
            improvements.append(best_improvement)
            neighbors.append(neighbor_ids)
        return np.array(improvements),relax_atom_id,relax_to_id,neighbors
    
    def update_cn_improvements(self,atoms,ids,atom_ids,free_positions,CNs,improvements,relax_to_id,relax_atom_id,neighbors):

        # Loop through free positions
        for i in ids:

            dists = np.linalg.norm(atoms.positions - free_positions[i],axis=1)

            nb_mask = dists<self.nb_dist

            n_nb = np.sum(nb_mask)
            if n_nb==0:
                relax_atom_id[i] = None
                relax_to_id[i] = i
                improvements[i] = -100
                neighbors[i] = []
                continue


            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            # Change in CN
            d_CN = (n_nb - 1) - CNs[neighbor_ids]

            # Get improvement from moving a neighbor into the free position
            best_improvement = np.max(d_CN)



            # Move a neighbor into free position if it improves CN
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_CN)==best_improvement])

            # Save the neighbor to move, to where and its improvement
            relax_atom_id[i]=idx
            relax_to_id[i] = i
            improvements[i] = best_improvement
            neighbors[i] = neighbor_ids
        


        return improvements, relax_atom_id,relax_to_id,neighbors


    def relax_particle_batch_cn(self,atoms,cn,nl, free_positions,trajectory):
        np.random.shuffle(free_positions)
        atom_ids=np.arange(len(atoms))

        improvements,relax_atom_id,relax_to_id,neighbors = self.get_cn_improvements(atoms,free_positions,cn,atom_ids)
        
        sort_ids = np.argsort(-improvements)
        relax_mask = improvements[sort_ids]>0
        relax_atom_id = np.array(relax_atom_id,dtype='object')
        relax_to_id = np.array(relax_to_id,dtype='object')

        while np.any(relax_mask):


            atoms_to_relax = []
            fps_to_occupy = []
            update_atom_ids = []
            
            for i in sort_ids[relax_mask]:
                
                if np.any(np.isin(neighbors[i],update_atom_ids)):
                    pass
                else:
                    atoms_to_relax.append(relax_atom_id[i])
                    fps_to_occupy.append(relax_to_id[i])
                    update_atom_ids = np.unique(np.concatenate((update_atom_ids,neighbors[i],nl.get_neighbors(relax_atom_id[i])[0])))


            atoms_to_relax = np.array(atoms_to_relax,dtype=int)
            fps_to_occupy = np.array(fps_to_occupy,dtype=int)


            new_free_positions = atoms.positions[atoms_to_relax].copy()

            atoms.positions[atoms_to_relax] = free_positions[fps_to_occupy].copy()

            free_positions[fps_to_occupy] = new_free_positions.copy()

            update_atom_ids = np.asarray(update_atom_ids,dtype=int)
            cn,nl = self.update_cn(atoms,update_atom_ids,nl,cn)
            
            if trajectory is not None:
                trajectory.write(atoms)

            affected_fp_nfp = np.concatenate([np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0] for new_free_position in new_free_positions])
            affected_fp_ofp = np.concatenate([np.where(np.linalg.norm(free_positions-atoms.positions[atom_to_relax],axis=1)<=self.nb_dist)[0] for atom_to_relax in atoms_to_relax])
            
            update_fp_ids = np.unique(np.append(affected_fp_nfp,affected_fp_ofp))

            improvements, relax_atom_id,relax_to_id,neighbors = self.update_cn_improvements(atoms,update_fp_ids,atom_ids,free_positions,cn,improvements,relax_to_id,relax_atom_id,neighbors)
            

            
            sort_ids = np.argsort(-improvements)
            relax_mask = improvements[sort_ids]>0            
            

        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]

        return atoms, cn, nl, free_positions
    

    def relax_particle_single_cn(self,atoms,cn,nl, free_positions,trajectory):
        np.random.shuffle(free_positions)
        atom_ids=np.arange(len(atoms))
        
        improvements,relax_atom_id,relax_to_id,neighbors = self.get_cn_improvements(atoms,free_positions,cn,atom_ids)

        sort_ids = np.argsort(-improvements)

        while np.any(np.round(improvements,decimals=12)>0):
            
            best_idx = sort_ids[0]

            old_neighbors = nl.get_neighbors(relax_atom_id[best_idx])[0]

            new_free_position = atoms.positions[relax_atom_id[best_idx]].copy()

            atoms.positions[relax_atom_id[best_idx]] = free_positions[relax_to_id[best_idx]].copy()

            free_positions[relax_to_id[best_idx]] = new_free_position.copy()
            
            update_ids = np.unique(np.append(old_neighbors,neighbors[best_idx]))

            cn,nl = self.update_cn(atoms,update_ids,nl,cn)
            
            
            if trajectory is not None:
                trajectory.write(atoms)


            update_ids = np.unique(np.append(np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0], np.where(np.linalg.norm(free_positions-atoms.positions[relax_atom_id[best_idx]],axis=1)<=self.nb_dist)[0]))
            
            improvements, relax_atom_id,relax_to_id,neighbors = self.update_cn_improvements(atoms,update_ids,atom_ids,free_positions,cn,improvements,relax_to_id,relax_atom_id,neighbors)
            sort_ids = np.argsort(-improvements)

            
        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]
        

        return atoms, cn, nl, free_positions



    def get_dE(self,symbols,features):
        dE = np.zeros(len(features))
        for metal in self.metals:
            mask = symbols==metal
            dE[mask] = self.regressor.predict(metal,features[mask])
        return dE
    
    def get_Udiss(self,symbols,dE):
        # initiate array
        U_diss = np.zeros_like(dE)

        # Calculate U_diss by each target metal
        for metal in self.metals:
            mask = symbols==metal
            U_diss[mask] = np.min(dE[mask]/self.U_dict[metal]['n'] + self.U_dict[metal]['U'],axis=0)
        return U_diss


    def dissolve_atoms_batch(self,applied_potential,atoms=None,relax_func=None,traj_file=None,return_trajectory=False):

        assert self.particle is not None or self.atoms is not None, 'Must provide atoms object or make particle using Dissolver.make_particle()'

        assert self.regressor is not None, 'A regressor must be set'

        assert self.feature_generator is not None, 'A feature generator must be set'

        if atoms is None:
            atoms = self.particle.copy()

        free_positions=np.empty((0,3))
        removed_atoms = []

        if traj_file is not None:
            trajectory = Trajectory(traj_file,'w')
            trajectory.write(atoms)
        else:
            trajectory=None
     
        if return_trajectory:
            trajectory_list = [atoms.copy()]


        # Symbol and IDs of atoms
        symbols = np.array(atoms.get_chemical_symbols())
        ids = np.arange(len(atoms))

        # Simulation
        # Get coordination numbers of all atoms
        coordination_numbers,nl = self.get_coordination_numbers(atoms,return_neighborlist=True)

        # Mask of atoms on surface to calculate
        surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
        surface_ids = ids[surface_mask]
        surface_symbols = np.array(symbols)[surface_mask]

        # get featuers of surface atoms
        features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

        # dEs[surface_mask] = self.get_dE(surface_symbols,features)
        dEs = self.get_dE(surface_symbols,features)

        while len(atoms)>0:
            
            if len(dEs)>0:
                # get dissolution potentials
                dissolution_potentials = self.get_Udiss(surface_symbols,dEs)
            
                # Dissolve atoms where U_diss<U
                remove_ids = surface_ids[dissolution_potentials<applied_potential]
            else:
                remove_ids = np.empty(0,dtype=int)

            # Dissolve atoms with lower CN than included in model
            if np.any(coordination_numbers<self.min_cn):
                remove_ids = np.append(remove_ids,ids[coordination_numbers<3])

            if len(remove_ids)==0:
                # Break simulation if there are no atoms to dissolve
                break
            else:
                # List the atoms to dissolve's position as free
                free_positions= np.vstack((free_positions,np.atleast_2d(atoms.positions[remove_ids])))
                # Save symbols of atoms to dissolve
                removed_atoms.append(symbols[remove_ids])
                # Dissolve atoms by deleting them from ASE atoms object
                del atoms[remove_ids]
                symbols = np.delete(symbols,remove_ids)
                ids = np.arange(len(atoms))

                if traj_file is not None:
                    trajectory.write(atoms)

                # Get new CN
                coordination_numbers,nl = self.get_coordination_numbers(atoms,True)


                # relax particle?
                if relax_func is not None:
                    atoms, coordination_numbers, nl, free_positions  = relax_func(atoms,coordination_numbers,nl,free_positions,trajectory)


                # Mask of atoms on surface to calculate
                surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                # dEs = np.zeros(len(atoms)) + np.inf
                
                
                surface_ids = ids[surface_mask]
                surface_symbols = symbols[surface_mask]

                if np.sum(surface_mask)>0:
                    # get featuers of surface atoms
                    features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

                    # dEs[surface_mask] = self.get_dE(surface_symbols,features)
                    dEs = self.get_dE(surface_symbols,features)
                else:
                    dEs = []
                
                if return_trajectory:
                    trajectory_list.append(atoms.copy())
                    

        if return_trajectory:
            return atoms,removed_atoms,trajectory_list
        else:
            if len(removed_atoms)>0:
                removed_atoms = np.concatenate(removed_atoms)
            return atoms, removed_atoms
    

    def dissolve_atoms_single(self,applied_potential,atoms=None,relax_func=None,traj_file=None,return_trajectory=False):

        assert self.particle is not None or self.atoms is not None, 'Must provide atoms object or make particle using Dissolver.make_particle()'

        assert self.regressor is not None, 'A regressor must be set'

        assert self.feature_generator is not None, 'A feature generator must be set'

        if atoms is None:
            atoms = self.particle.copy()

        free_positions=np.empty((0,3))
        removed_atoms = np.empty(0,dtype=str)

        if traj_file is not None:
            trajectory = Trajectory(traj_file,'w')
            trajectory.write(atoms)
        else:
            trajectory=None
     
        if return_trajectory:
            trajectory_list = [atoms.copy()]


        # Symbol and IDs of atoms
        symbols = np.array(atoms.get_chemical_symbols())
        ids = np.arange(len(atoms))

        # Simulation
        # Get coordination numbers of all atoms
        coordination_numbers,nl = self.get_coordination_numbers(atoms,return_neighborlist=True)

        # Mask of atoms on surface to calculate
        surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
        surface_ids = ids[surface_mask]
        surface_symbols = np.array(symbols)[surface_mask]


        # get featuers of surface atoms
        features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

        dEs = self.get_dE(surface_symbols,features)

        while len(atoms)>0:
            
            if len(dEs)>0:
                # get dissolution potentials
                dissolution_potentials = self.get_Udiss(surface_symbols,dEs)

                min_Udiss = np.min(dissolution_potentials)

                if min_Udiss>=applied_potential:
                    break

                # Dissolve atoms where U_diss<U
                remove_ids = [np.random.choice(surface_ids[dissolution_potentials==min_Udiss])]
            else:
                remove_ids = np.empty(0,dtype=int)

            # Dissolve atoms with lower CN than included in model
            if np.any(coordination_numbers<self.min_cn):
                remove_ids = np.astype(np.append(remove_ids,ids[coordination_numbers<3]),int)
            if len(remove_ids)==0:
                # Break simulation if there are no atoms to dissolve
                break
            else:
                # List the atoms to dissolve's position as free
                if relax_func is not None:
                    free_positions= np.vstack((free_positions,np.atleast_2d(atoms.positions[remove_ids])))
                # Save symbols of atoms to dissolve
                removed_atoms = np.append(removed_atoms,symbols[remove_ids])
                # Dissolve atoms by deleting them from ASE atoms object
                del atoms[remove_ids]

                symbols = np.delete(symbols,remove_ids)
                ids = np.arange(len(atoms))

                if traj_file is not None:
                    trajectory.write(atoms)

                # Get new CN
                coordination_numbers,nl = self.get_coordination_numbers(atoms,True)


                # relax particle?
                if relax_func is not None:
                    atoms, coordination_numbers, nl, free_positions = relax_func(atoms,coordination_numbers,nl,free_positions,trajectory)

                
                # Mask of atoms on surface to calculate
                surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                
                
                surface_ids = ids[surface_mask]
                surface_symbols = symbols[surface_mask]

                if np.any(surface_mask):
                    # get featuers of surface atoms
                    features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

                    dEs = self.get_dE(surface_symbols,features)
                else:
                    dEs = []

                
                if return_trajectory:
                    trajectory_list.append(atoms.copy())
                    

        if return_trajectory:
            return atoms,removed_atoms,trajectory_list
        else:
            return atoms, removed_atoms