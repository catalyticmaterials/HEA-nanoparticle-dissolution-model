import numpy as np

def get_site_pos_in_xy(pos_1st):
    # Rearange positions so that they are ordered by layer
    pos_1st_ = np.array([pos_1st[(9*3*i)+j*3:3+3*9*i+j*3] for j in range(9) for i in range(3)])
    pos_1st_ = pos_1st_.reshape(81,3)
    #print(pos_1st_)
    grid = pos_1st_.reshape(9,9,-1)
    #print(grid)
    grid = np.pad(grid,pad_width=((1,1),(1,1),(0,0)),mode="wrap")
    
    fcc_sites =  (grid[1:-1,1:-1] + grid[1:-1,2:] + grid[2:,1:-1])/3
    
    hcp_sites = (grid[1:-1,1:-1] + grid[2:,:-2] +grid[2:,1:-1])/3
    
    bridge_sites1 = (grid[1:-1,1:-1] + grid[1:-1,2:])/2
    bridge_sites2 = (grid[1:-1,1:-1] + grid[2:,1:-1])/2
    bridge_sites3 = (grid[1:-1,1:-1] + grid[2:,:-2])/2
    
    
    # cut off ends as we are only interested in sites in the middle 3x3 atoms anyway
    fcc_sites = fcc_sites[2:-2,2:-2]
    hcp_sites = hcp_sites[2:-2,2:-2]
    bridge_sites1 = bridge_sites1[2:-2,2:-2]
    bridge_sites2 = bridge_sites2[2:-2,2:-2]
    bridge_sites3 = bridge_sites3[2:-2,2:-2]
    ontop_sites = np.copy(grid[3:-3,3:-3]).reshape(-1,3)
    
    bridge_sites = np.vstack([bridge_sites1.reshape(-1,3),bridge_sites2.reshape(-1,3),bridge_sites3.reshape(-1,3)])
    return fcc_sites.reshape(-1,3), hcp_sites.reshape(-1,3), bridge_sites,ontop_sites

def get_nearest_sites_in_xy(fcc,hcp,bridge,ontop,ads):
    fcc_dist = np.sum((fcc[:,:2]-ads[:2])**2,axis=1)
    hcp_dist = np.sum((hcp[:,:2]-ads[:2])**2,axis=1)
    bridge_dist = np.sum((bridge[:,:2]-ads[:2])**2,axis=1)
    ontop_dist = np.sum((ontop[:,:2]-ads[:2])**2,axis=1)
    
    min_ids = [np.argmin(dist) for dist in (fcc_dist,hcp_dist,bridge_dist,ontop_dist)]
    min_dists = [dist[min_ids[i]] for i,dist in enumerate((fcc_dist,hcp_dist,bridge_dist,ontop_dist))]
    
    site_str = ["fcc","hcp","bridge","ontop"]
    
    nearest_site_type = site_str[np.argmin(min_dists)]
    
    return nearest_site_type, min_ids



def check_site(atoms,ads,n_atoms_site=3,return_positions=False):
    # Repeat atoms object
    atoms_3x3 = atoms.repeat((3,3,1))
    
    # Get chemical symbols of atoms
    symbols = np.asarray(atoms_3x3.get_chemical_symbols())
    
    # Get index of central adsorbate atom
    idx_ads= np.nonzero(symbols == ads)[0][4]
    
    # Get position of hydrogen atom
    pos_ads = atoms_3x3.positions[idx_ads]
    
    # Get indices of the 1st layer atoms
    ids_1st = np.array([atom.index for atom in atoms_3x3 if atom.tag == 1])
    
    # Get posittion of atoms in 1st layer
    pos_1st=atoms_3x3.positions[ids_1st]
    
    #Get position of each site type in xy-plane from the 1st layer atom position
    sites_xy = get_site_pos_in_xy(pos_1st)
    # fcc_sites_pos, hcp_sites_pos, bridge_sites_pos, ontop_sites_pos = sites_xy
    #Get site type and min dist site id.
    nearest_site, min_dist_ids = get_nearest_sites_in_xy(*sites_xy,pos_ads)

    site_ids = ids_1st[np.argsort(np.linalg.norm(pos_ads - pos_1st,axis=1))[:n_atoms_site]]

    site_symbols = symbols[site_ids]

    if return_positions:
        return sites_xy, pos_ads, min_dist_ids ,nearest_site, site_ids, site_symbols

    return nearest_site, site_ids, site_symbols