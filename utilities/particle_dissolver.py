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

    def get_improvements(self,atoms,free_positions,dEs,atom_ids,symbols,cn):
        # Get best improvement for each unoccupied position
        # discard_free_pos = []
        relax_to_id = []
        relax_atom_id = []
        improvements = []
        neighbors = []
        # Loop through free positions
        for i,free_position in enumerate(free_positions):
            
            dists = np.linalg.norm(atoms.positions - free_position,axis=1)

            nb_mask = dists<self.nb_dist

            n_nb = np.sum(nb_mask)

            

            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            if n_nb<(self.min_cn+1):
                relax_atom_id.append(np.nan)
                relax_to_id.append(i)
                improvements.append(np.inf)
                neighbors.append(neighbor_ids)
                continue
            elif n_nb>(self.max_cn+1):
                relax_atom_id.append(np.random.choice(neighbor_ids[dEs[neighbor_ids]<np.inf]))
                relax_to_id.append(i)
                improvements.append(-np.inf)
                neighbors.append(neighbor_ids)
                continue

            

            # neighbor_symbols = symbols[neighbor_ids]

            # neighbor_feature = [np.sum(neighbor_symbols==metal) for metal in self.metals]
            # cn_feature = np.zeros(7)
            # cn_feature[n_nb-4] = 1
            d_dE = []
            for neighbor in neighbor_ids:
                # neighbor=13
                # neighbor_feature_copy = deepcopy(neighbor_feature)
                target_metal = symbols[neighbor]
                # neighbor_feature_copy[metals.index(target_metal)] -=1

                # feature = np.append(cn_feature, neighbor_feature_copy/(n_nb-1))
                
                
                # d_dE.append(dEs[neighbor] - self.regressor.predict(target_metal,feature))

                old_cn_idx = cn[neighbor]-3
                if old_cn_idx>6:
                    d_dE.append(np.inf)
                elif old_cn_idx <0:
                    d_dE.append(-np.inf)
                else:
                    new_cn_idx = n_nb-1 -3
                    d_dE.append(self.regressor.parameters[target_metal][old_cn_idx] - self.regressor.parameters[target_metal][new_cn_idx])

                

            # Get improvement from moving a neighbor into the free position
            best_improvement = np.min(d_dE)

            # Move a neighbor into free position if it improves CN
            # if best_improvement<0:
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_dE)==best_improvement])
    
            # Save the neighbor to move, to where and its improvement
            relax_atom_id.append(idx)
            relax_to_id.append(i)
            improvements.append(best_improvement)
            neighbors.append(neighbor_ids)
        return np.array(improvements),relax_atom_id,relax_to_id,neighbors
    
    def update_improvements(self,atoms,ids,atom_ids,symbols,free_positions,dEs,improvements,relax_to_id,relax_atom_id,neighbors,cn):
        
        # Get best improvement for each unoccupied position
        # discard_free_pos = []

        # Loop through free positions
        for i in ids:
            
            dists = np.linalg.norm(atoms.positions - free_positions[i],axis=1)

            nb_mask = dists<self.nb_dist

            n_nb = np.sum(nb_mask)

            # assert n_nb <=10, f'{n_nb} neigbors. Too high coordinated free position'


            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            if n_nb<(self.min_cn+1):
                relax_atom_id[i] = np.nan
                relax_to_id[i] = i
                improvements[i] = np.inf
                neighbors[i] = neighbor_ids
                continue
            elif n_nb>(self.max_cn+1):
                relax_atom_id[i] = np.random.choice(neighbor_ids[dEs[neighbor_ids]<np.inf])
                relax_to_id[i] = i
                improvements[i] = -np.inf
                neighbors[i] = neighbor_ids
                continue

            # neighbor_symbols = symbols[neighbor_ids]

            # neighbor_feature = [np.sum(neighbor_symbols==metal) for metal in self.metals]
            # cn_feature = np.zeros(7)
            # cn_feature[n_nb-self.min_cn-1] = 1
            new_cn_idx = n_nb-1 - 3
            d_dE = []
            for neighbor in neighbor_ids:
                # neighbor=13
                # neighbor_feature_copy = deepcopy(neighbor_feature)
                target_metal = symbols[neighbor]
                # neighbor_feature_copy[metals.index(target_metal)] -=1

                # feature = np.append(cn_feature, neighbor_feature_copy/(n_nb-1))
                
                # d_dE.append(dEs[neighbor] - self.regressor.predict(target_metal,feature))
                old_cn_idx = cn[neighbor]-3
                if old_cn_idx>6:
                    d_dE.append(np.inf)
                elif old_cn_idx <0:
                    d_dE.append(-np.inf)
                else:
                    
                    d_dE.append(self.regressor.parameters[target_metal][old_cn_idx] - self.regressor.parameters[target_metal][new_cn_idx])
                    
   
            # Get improvement from moving a neighbor into the free position
            best_improvement = np.min(d_dE)

            # Move a neighbor into free position if it improves CN
            # if best_improvement<0:
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_dE)==best_improvement])
    
            # Save the neighbor to move, to where and its improvement
            relax_atom_id[i]=idx
            relax_to_id[i] = i
            improvements[i] = best_improvement
            neighbors[i] = neighbor_ids
        


        return improvements, relax_atom_id,relax_to_id,neighbors

    def relax_particle_single(self,atoms,cn,nl,dEs,symbols, free_positions,trajectory):

        atom_ids=np.arange(len(atoms))
        
        improvements,relax_atom_id,relax_to_id,neighbors = self.get_improvements(atoms,free_positions,dEs,atom_ids,symbols,cn)
        # latest_moved_atoms = deque(maxlen=8)
        # moved_atoms = []
        sort_ids = np.argsort(improvements)
        # not_moved_mask = np.isin(np.array(relax_atom_id)[sort_ids],moved_atoms,invert=True)
        # while np.any(np.round(improvements[sort_ids][not_moved_mask],decimals=12)<0):
        while np.any(np.round(improvements[sort_ids],decimals=12)<0):
            
            # best_idx = sort_ids[not_moved_mask][0]
            best_idx = sort_ids[0]

            old_neighbors = nl.get_neighbors(relax_atom_id[best_idx])[0]

            new_free_position = atoms.positions[relax_atom_id[best_idx]].copy()

            atoms.positions[relax_atom_id[best_idx]] = free_positions[relax_to_id[best_idx]].copy()

            free_positions[relax_to_id[best_idx]] = new_free_position.copy()
            
            # moved_atoms.append(relax_atom_id[best_idx])
            
            update_ids = np.append(old_neighbors,neighbors[best_idx])

            cn,nl = self.update_cn(atoms,update_ids,nl,cn)

            # dEs = self.update_dE(dEs,update_ids,cn,nl,symbols)
            
            if trajectory is not None:
                trajectory.write(atoms)

            # latest_moved_atoms.append(relax_atom_id[best_idx])
            # if np.all(np.unique(latest_moved_atoms,return_counts=True)[1]==4):
            #     # print(latest_moved_atoms,np.all(np.unique(latest_moved_atoms,return_counts=True)[1]==4))
            #     # stop
            #     break

            update_ids = np.unique(np.append(np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0], np.where(np.linalg.norm(free_positions-atoms.positions[relax_atom_id[best_idx]],axis=1)<=self.nb_dist)[0]))
            
            improvements, relax_atom_id,relax_to_id,neighbors = self.update_improvements(atoms,update_ids,atom_ids,symbols,free_positions,dEs,improvements,relax_to_id,relax_atom_id,neighbors,cn)
            sort_ids = np.argsort(improvements)
            # not_moved_mask = np.isin(np.array(relax_atom_id)[sort_ids],moved_atoms,invert=True)
            

            
        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]
        

        return atoms, cn, nl, dEs, free_positions, True


    def relax_particle_batch(self,atoms,cn,nl,dEs,symbols, free_positions,trajectory):
        atom_ids=np.arange(len(atoms))
        improvements,relax_atom_id,relax_to_id,neighbors = self.get_improvements(atoms,free_positions,dEs,atom_ids,symbols,cn)

        sort_ids = np.argsort(improvements)

        relax_mask = np.round(improvements[sort_ids],decimals=12)<0
        relax_atom_id = np.array(relax_atom_id,dtype='object')
        relax_to_id = np.array(relax_to_id,dtype='object')
        # moved_atoms = np.empty(0)
        while np.any(relax_mask):
            
            atoms_to_relax = []
            fps_to_occupy = []
            # new_neighbors = []
            update_ids = []
            for i in sort_ids[relax_mask]:
                if np.any(np.isin(neighbors[i],update_ids)):
                    pass
                else:
                    atoms_to_relax.append(relax_atom_id[i])
                    fps_to_occupy.append(relax_to_id[i])
                    update_ids = np.unique(np.concatenate((update_ids,neighbors[i],nl.get_neighbors(relax_atom_id[i])[0])))


            # print(atoms_to_relax,improvements[sort_ids][relax_mask])
            atoms_to_relax = np.array(atoms_to_relax,dtype=int)

            fps_to_occupy = np.array(fps_to_occupy,dtype=int)
            # new_neighbors = np.concatenate(new_neighbors)

            # old_neighbors = np.concatenate([nl.get_neighbors(idx)[0] for idx in atoms_to_relax])

            new_free_positions = atoms.positions[atoms_to_relax].copy()

            atoms.positions[atoms_to_relax] = free_positions[fps_to_occupy].copy()

            free_positions[fps_to_occupy] = new_free_positions.copy()

            # moved_atoms = np.append(moved_atoms,atoms_to_relax)

            # update_ids = np.unique(np.append(old_neighbors,np.concatenate([neighbor_ids for relax_bool,neighbor_ids in zip(relax_mask,neighbors) if relax_bool])))
            # update_ids = np.unique(np.append(old_neighbors,new_neighbors))
            update_ids = np.asarray(update_ids,dtype=int)
            # old_cn = cn.copy()
            cn,nl = self.update_cn(atoms,update_ids,nl,cn)
            # cn,nl = self.get_coordination_numbers(atoms,True)
            # cn_diff = cn[atoms_to_relax]-old_cn[atoms_to_relax]
            # assert_mask =cn_diff>=0

            

            # dEs = self.update_dE(dEs,update_ids,cn,nl,symbols)
            
            if trajectory is not None:
                trajectory.write(atoms)
            # assert np.all(assert_mask), f'{atoms_to_relax[np.invert(assert_mask)]},{cn_diff[np.invert(assert_mask)]}'
            
            affected_fp_nfp = np.concatenate([np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0] for new_free_position in new_free_positions])
            affected_fp_ofp = np.concatenate([np.where(np.linalg.norm(free_positions-atoms.positions[atom_to_relax],axis=1)<=self.nb_dist)[0] for atom_to_relax in atoms_to_relax])
            
            update_ids = np.unique(np.append(affected_fp_nfp,affected_fp_ofp))

            improvements, relax_atom_id,relax_to_id,neighbors = self.update_improvements(atoms,update_ids,atom_ids,symbols,free_positions,dEs,improvements,relax_to_id,relax_atom_id,neighbors,cn)
            sort_ids = np.argsort(improvements)
            relax_mask = np.round(improvements[sort_ids],decimals=12)<0

            # moved_mask = np.isin(relax_atom_id[sort_ids],moved_atoms)
            # relax_mask[moved_mask] = False
            

            
        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]
        return atoms, cn, nl, dEs, free_positions, True
    

    def get_cn_improvements(self,atoms,free_positions,CNs,atom_ids,symbols):
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
                improvements.append(-np.inf)
                neighbors.append([])
                continue


            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            
            d_CN = (n_nb - 1) - CNs[neighbor_ids]
           
            # Get improvement from moving a neighbor into the free position
            best_improvement = np.max(d_CN)

            # Move a neighbor into free position if it improves CN
            # if best_improvement<0:
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_CN)==best_improvement])
    
            # Save the neighbor to move, to where and its improvement
            relax_atom_id.append(idx)
            relax_to_id.append(i)
            improvements.append(best_improvement)
            neighbors.append(neighbor_ids)
        return np.array(improvements),relax_atom_id,relax_to_id,neighbors
    
    def update_cn_improvements(self,atoms,ids,atom_ids,free_positions,CNs,improvements,relax_to_id,relax_atom_id,neighbors):
        
        # Get best improvement for each unoccupied position
        # discard_free_pos = []

        # Loop through free positions
        for i in ids:

            dists = np.linalg.norm(atoms.positions - free_positions[i],axis=1)

            nb_mask = dists<self.nb_dist

            n_nb = np.sum(nb_mask)
            if n_nb==0:
                relax_atom_id[i] = None
                relax_to_id[i] = i
                improvements[i] -np.inf
                neighbors[i] = []
                continue


            # Get neighbors to free position 
            neighbor_ids = atom_ids[nb_mask]

            
            d_CN = (n_nb - 1) - CNs[neighbor_ids]

            # Get improvement from moving a neighbor into the free position

            best_improvement = np.max(d_CN)



            # Move a neighbor into free position if it improves CN
            # if best_improvement<0:
            # Choose neighbor to move by max improvement. Choose randomly out of all with max improvement.
            idx = np.random.choice(neighbor_ids[np.array(d_CN)==best_improvement])
    
            # Save the neighbor to move, to where and its improvement
            relax_atom_id[i]=idx
            relax_to_id[i] = i
            improvements[i] = best_improvement
            neighbors[i] = neighbor_ids
        


        return improvements, relax_atom_id,relax_to_id,neighbors


    def relax_particle_batch_cn(self,atoms,cn,nl,dEs,symbols, free_positions,trajectory):
        atom_ids=np.arange(len(atoms))

        improvements,relax_atom_id,relax_to_id,neighbors = self.get_cn_improvements(atoms,free_positions,cn,atom_ids,symbols)
        
        sort_ids = np.argsort(-improvements)
        relax_mask = improvements[sort_ids]>0
        relax_atom_id = np.array(relax_atom_id,dtype='object')
        relax_to_id = np.array(relax_to_id,dtype='object')
        # moved_atoms = np.empty(0)
        while np.any(relax_mask):


            atoms_to_relax = []
            fps_to_occupy = []
            # new_neighbors = []
            update_ids = []
            for i in sort_ids[relax_mask]:

                if np.any(np.isin(neighbors[i],update_ids)):
                    pass
                else:
                    atoms_to_relax.append(relax_atom_id[i])
                    fps_to_occupy.append(relax_to_id[i])
                    update_ids = np.unique(np.concatenate((update_ids,neighbors[i],nl.get_neighbors(relax_atom_id[i])[0])))


            # print(atoms_to_relax,improvements[sort_ids][relax_mask])
            atoms_to_relax = np.array(atoms_to_relax,dtype=int)
            fps_to_occupy = np.array(fps_to_occupy,dtype=int)
            # new_neighbors = np.concatenate(new_neighbors)

            # old_neighbors = np.concatenate([nl.get_neighbors(idx)[0] for idx in atoms_to_relax])

            new_free_positions = atoms.positions[atoms_to_relax].copy()

            atoms.positions[atoms_to_relax] = free_positions[fps_to_occupy].copy()

            free_positions[fps_to_occupy] = new_free_positions.copy()

            # moved_atoms = np.append(moved_atoms,atoms_to_relax)

            # update_ids = np.unique(np.append(old_neighbors,np.concatenate([neighbor_ids for relax_bool,neighbor_ids in zip(relax_mask,neighbors) if relax_bool])))
            # update_ids = np.unique(np.append(old_neighbors,new_neighbors))
            update_ids = np.asarray(update_ids,dtype=int)
            cn,nl = self.update_cn(atoms,update_ids,nl,cn)
            
            if trajectory is not None:
                trajectory.write(atoms)

            affected_fp_nfp = np.concatenate([np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0] for new_free_position in new_free_positions])
            affected_fp_ofp = np.concatenate([np.where(np.linalg.norm(free_positions-atoms.positions[atom_to_relax],axis=1)<=self.nb_dist)[0] for atom_to_relax in atoms_to_relax])
            
            update_ids = np.unique(np.append(affected_fp_nfp,affected_fp_ofp))

            improvements, relax_atom_id,relax_to_id,neighbors = self.update_cn_improvements(atoms,update_ids,atom_ids,free_positions,cn,improvements,relax_to_id,relax_atom_id,neighbors)
            
            sort_ids = np.argsort(-improvements)
            relax_mask = improvements[sort_ids]>0
            # moved_mask = np.isin(relax_atom_id[sort_ids],moved_atoms)
            # relax_mask[moved_mask] = False
            

        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]

        return atoms, cn, nl, dEs, free_positions, True
    

    def relax_particle_single_cn(self,atoms,cn,nl,dEs,symbols, free_positions,trajectory):

        atom_ids=np.arange(len(atoms))
        
        improvements,relax_atom_id,relax_to_id,neighbors = self.get_cn_improvements(atoms,free_positions,cn,atom_ids,symbols)
        # latest_moved_atoms = deque(maxlen=8)
        # moved_atoms = []
        sort_ids = np.argsort(-improvements)
        # not_moved_mask = np.isin(np.array(relax_atom_id)[sort_ids],moved_atoms,invert=True)
        # while np.any(np.round(improvements[sort_ids][not_moved_mask],decimals=12)<0):
        while np.any(np.round(improvements,decimals=12)>0):
            
            # best_idx = sort_ids[not_moved_mask][0]
            best_idx = sort_ids[0]

            old_neighbors = nl.get_neighbors(relax_atom_id[best_idx])[0]

            new_free_position = atoms.positions[relax_atom_id[best_idx]].copy()

            atoms.positions[relax_atom_id[best_idx]] = free_positions[relax_to_id[best_idx]].copy()

            free_positions[relax_to_id[best_idx]] = new_free_position.copy()
            
            # moved_atoms.append(relax_atom_id[best_idx])
            
            update_ids = np.unique(np.append(old_neighbors,neighbors[best_idx]))

            cn,nl = self.update_cn(atoms,update_ids,nl,cn)
            
            # dEs = self.update_dE(dEs,update_ids,cn,nl,symbols)
            
            if trajectory is not None:
                trajectory.write(atoms)

            # latest_moved_atoms.append(relax_atom_id[best_idx])
            # if np.all(np.unique(latest_moved_atoms,return_counts=True)[1]==4):
            #     # print(latest_moved_atoms,np.all(np.unique(latest_moved_atoms,return_counts=True)[1]==4))
            #     # stop
            #     break

            update_ids = np.unique(np.append(np.where(np.linalg.norm(free_positions-new_free_position,axis=1)<=self.nb_dist)[0], np.where(np.linalg.norm(free_positions-atoms.positions[relax_atom_id[best_idx]],axis=1)<=self.nb_dist)[0]))
            
            improvements, relax_atom_id,relax_to_id,neighbors = self.update_cn_improvements(atoms,update_ids,atom_ids,free_positions,cn,improvements,relax_to_id,relax_atom_id,neighbors)
            sort_ids = np.argsort(-improvements)
            # not_moved_mask = np.isin(np.array(relax_atom_id)[sort_ids],moved_atoms,invert=True)
            

            
        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]
        

        return atoms, cn, nl, dEs, free_positions, True



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

        dEs = np.zeros(len(atoms)) + np.inf

        # get featuers of surface atoms
        features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

        dEs[surface_mask] = self.get_dE(surface_symbols,features)

        while len(atoms)>0:
            
            # get dissolution potentials
            dissolution_potentials = self.get_Udiss(surface_symbols,dEs[surface_mask])
        
            # Dissolve atoms where U_diss<U
            remove_ids = surface_ids[dissolution_potentials<applied_potential]

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
                removed_atoms = np.append(removed_atoms,symbols[remove_ids])
                # Dissolve atoms by deleting them from ASE atoms object
                del atoms[remove_ids]

                symbols = np.delete(symbols,remove_ids)
                ids = np.arange(len(atoms))

                if traj_file is not None:
                    trajectory.write(atoms)


                coordination_numbers,nl = self.get_coordination_numbers(atoms,True)

                # Mask of atoms on surface to calculate
                surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                dEs = np.zeros(len(atoms)) + np.inf
                
                
                surface_ids = ids[surface_mask]
                surface_symbols = symbols[surface_mask]

                if np.sum(surface_mask)>0:
                    # get featuers of surface atoms
                    features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

                    dEs[surface_mask] = self.get_dE(surface_symbols,features)


                # relax particle?
                if relax_func is not None:
                    np.random.shuffle(free_positions)
                    atoms, coordination_numbers, nl, dEs, free_positions, update_dE = relax_func(atoms,coordination_numbers,nl,dEs,symbols,free_positions,trajectory)
                    # surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                    # surface_ids = ids[surface_mask]
                    # surface_symbols = symbols[surface_mask]
                    # if update_dE:
                    #     # get featuers of surface atoms
                    #     features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)
                    #     dEs[surface_mask] = self.get_dE(surface_symbols,features)

                # Mask of atoms on surface to calculate
                surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                dEs = np.zeros(len(atoms)) + np.inf
                
                
                surface_ids = ids[surface_mask]
                surface_symbols = symbols[surface_mask]

                if np.sum(surface_mask)>0:
                    # get featuers of surface atoms
                    features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

                    dEs[surface_mask] = self.get_dE(surface_symbols,features)
                    

        if return_trajectory:
            # trajectory_list.append(atoms.copy())
            return atoms,removed_atoms,trajectory_list
        else:
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

        dEs = np.zeros(len(atoms)) + np.inf

        # get featuers of surface atoms
        features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

        dEs[surface_mask] = self.get_dE(surface_symbols,features)

        while len(atoms)>0:
            
            if np.any(surface_mask):
                # get dissolution potentials
                dissolution_potentials = self.get_Udiss(surface_symbols,dEs[surface_mask])

                min_Udiss = np.min(dissolution_potentials)

                if min_Udiss>=applied_potential:
                    break

                # Dissolve atoms where U_diss<U
                remove_ids = [np.random.choice(surface_ids[dissolution_potentials==min_Udiss])]
            else:
                remove_ids = []

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


                coordination_numbers,nl = self.get_coordination_numbers(atoms,True)


                # relax particle?
                if relax_func is not None:
                    np.random.shuffle(free_positions)
                    atoms, coordination_numbers, nl, dEs, free_positions, update_dE = relax_func(atoms,coordination_numbers,nl,dEs,symbols,free_positions,trajectory)
                    # surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                    # surface_ids = ids[surface_mask]
                    # surface_symbols = symbols[surface_mask]
                    # if update_dE:
                    #     # get featuers of surface atoms
                    #     features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)
                    #     dEs[surface_mask] = self.get_dE(surface_symbols,features)

                
                # Mask of atoms on surface to calculate
                surface_mask = (coordination_numbers<=self.max_cn)*(coordination_numbers>=self.min_cn)
                dEs = np.zeros(len(atoms)) + np.inf
                
                
                surface_ids = ids[surface_mask]
                surface_symbols = symbols[surface_mask]

                if np.any(surface_mask):
                    # get featuers of surface atoms
                    features = self.feature_generator(nl,surface_ids,symbols,coordination_numbers)

                    dEs[surface_mask] = self.get_dE(surface_symbols,features)
                    

        if return_trajectory:
            # trajectory_list.append(atoms.copy())
            return atoms,removed_atoms,trajectory_list
        else:
            return atoms, removed_atoms
        


    
    def get_total_dE(self,atoms_):
        atoms = atoms_.copy()
        tot = 0

        cn,nl = self.get_coordination_numbers(atoms,True)
        surface_mask = (cn<=9)*(cn>=3)
        symbols = np.array(atoms.get_chemical_symbols())
        ids = np.where(surface_mask)[0]
        features = self.feature_generator(nl,ids,symbols,cn)
        dEs = np.zeros(len(atoms)) + np.inf
        dEs[surface_mask] = self.get_dE(symbols[surface_mask],features)
        
        while True:
            
            

            tot+=np.min(dEs)

            remove_id = np.argmin(dEs)
        
            position = atoms.positions[remove_id]

            del atoms[remove_id]

            symbols = np.delete(symbols,remove_id)
            dEs = np.delete(dEs,remove_id)
            cn = np.delete(cn,remove_id)

            neighbors = np.where(np.linalg.norm(atoms.positions - position,axis=1)<=self.nb_dist)[0]

            nl = NeighborList([self.nb_dist/2]*len(atoms))
            nl.update(atoms)
            cn[neighbors] = [len(nl.get_neighbors(idx)[0]) for idx in neighbors]

            if np.any(cn<3):
                remove_ids = np.where(cn<3)[0]

                del atoms[remove_ids]

                cn,nl = self.get_coordination_numbers(atoms,True)
                surface_mask = (cn<=9)*(cn>=3)
                if np.all(surface_mask==False):
                    break
                symbols = np.array(atoms.get_chemical_symbols())
                ids = np.where(surface_mask)[0]
                features = self.feature_generator(nl,ids,symbols,cn)
                dEs = np.zeros(len(atoms)) + np.inf
                dEs[surface_mask] = self.get_dE(symbols[surface_mask],features)
            else:
                surface_mask = (cn<=9)*(cn>=3)

                if np.all(surface_mask==False):
                    break

                dEs = self.update_dE(dEs,neighbors,cn,nl,symbols)

 
        return tot
    


    def relax_mc(self,atoms,cn,nl,dEs,symbols, free_positions,trajectory):
        
        dE_total = self.get_total_dE(atoms)

        count = 0
        while count<20:
            free_position_id = np.random.choice(range(len(free_positions)))
            free_position = free_positions[free_position_id]

            nb_mask = np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist

            neighbors = np.where(nb_mask)[0]
            if len(neighbors)==0: continue
            neighbor = np.random.choice(neighbors)

            atoms_new = atoms.copy()

            atoms_new.positions[neighbor] = free_position


            dE_total_new = self.get_total_dE(atoms_new)

            dE_total_diff = dE_total_new - dE_total

            if dE_total_diff>0 or np.random.random()<=np.exp(dE_total_diff/(kB*T)):

                free_positions[free_position_id] = atoms.positions[neighbor].copy()

                atoms = atoms_new.copy()

                count = 0

                dE_total = dE_total_new

                if trajectory is not None:
                    trajectory.write(atoms)

                
            else:
                count+=1
                
            print(count,dE_total)
            
            
        mask = np.array([np.sum(np.linalg.norm(atoms.positions - free_position,axis=1)<=self.nb_dist)>=3 for free_position in free_positions])
        free_positions = free_positions[mask]

        return atoms, cn, nl, dEs, free_positions, True