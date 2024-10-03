from .Simplex import Simplex
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import itertools as it
import scipy
#import torch 
from copy import deepcopy

def count_elements(elements, n_elems):
	count = np.zeros(n_elems, dtype=int)
	for elem in elements:
		count[elem] += 1
	return count

def get_molar_fractions(step_size, n_elems, total=1., return_number_of_molar_fractions=False):
	'Get all molar fractions with the given step size'
	
	interval = int(total/step_size)
	n_combs = scipy.special.comb(n_elems+interval-1, interval, exact=True)
	
	if return_number_of_molar_fractions:
		return n_combs
		
	counts = np.zeros((n_combs, n_elems), dtype=int)

	for i, comb in enumerate(it.combinations_with_replacement(range(n_elems), interval)):
		counts[i] = count_elements(comb, n_elems)

	return counts*step_size

def get_random_molar_fractions(n_elems, n_molar_fractions=1, random_state=None):
	'Get ´size´ random molar fractions of ´n_elems´ elements'
	if random_state is not None:
		np.random.seed(random_state)

	fs = np.random.rand(n_molar_fractions, n_elems)
	return fs / np.sum(fs, axis=1)[:, None]

	
def get_composition(f, metals, return_latex=False, saferound=True):
	
	# Make into numpy and convert to atomic percent
	f = np.asarray(f)*100
	
	if saferound:
		# Round while maintaining the sum, the iteround module may need
		# to be installed manually from pypi: "pip3 install iteround"
		import iteround
		f = iteround.saferound(f, 0)
	
	if return_latex:
		# Return string in latex format with numbers as subscripts
		return ''.join(['$\\rm {0}_{{{1}}}$'.format(m,f0) for m,f0 in\
			zip(metals, map('{:.0f}'.format, f)) if float(f0) > 0.])
	else:
		# Return composition as plain text
		return ''.join([''.join([m, f0]) for m,f0 in\
			zip(metals, map('{:.0f}'.format, f)) if float(f0) > 0.])

def molar_fractions_to_cartesians(fs):
	
	# Make into numpy
	fs = np.asarray(fs)

	if fs.ndim == 1:
		fs = np.reshape(fs, (1, -1))

	# Get vertices of the multidimensional simplex
	n_elems = fs.shape[1]
	vertices = Simplex(n_elems).get_vertices()
	
	# Get cartesian coordinates corresponding to the molar fractions
	return fs @ vertices

def cartesians_to_molar_fractions(rs):
	
	# Make into numpy
	rs = np.asarray(rs)
	
	# Add column of ones to ´rs´ to account for the restriction
	# that the molar fractions must sum to unity
	rs = np.concatenate((rs, np.ones((rs.shape[0], 1))), axis=1)
	
	# Get vertices of the multidimensional simplex
	n_elems = rs.shape[1]
	vertices = Simplex(n_elems).get_vertices()
	
	# Add column of ones to ´vertices´ to account for the restriction
	# that the molar fractions must sum to unity
	vertices = np.concatenate((vertices, np.ones((vertices.shape[0], 1))), axis=1)
	
	# Get molar fractions corresponding to the cartesian coordinates
	# r = fV <=> r^T = (fV)^T = V^T f^T (i.e. on the form Ax = b that np.linalg.solve takes as input)
	return np.linalg.solve(vertices.T, rs.T).T

def get_molar_fractions_around(f, step_size, total=1., eps=1e-10):
	'Get all molar fractions with the given step size around the given molar fraction'	
	fs = []	
	n_elems = len(f)
	for pair, ids in zip(it.permutations(f, 2), it.permutations(range(n_elems), 2)):
	
		# Get molar fractions and their ids
		f0, f1 = pair
		id0, id1 = ids
		
		# Increment one molar fraction and decrement the other
		f0_new = f0 + (step_size - eps)
		f1_new = f1 - (step_size - eps)
		
		# Ignore if the new molar fractions are not between 0 and 1
		if f0_new <= total and f1_new >= 0.:
			
			# Make new molar fraction
			f_new = deepcopy(f)
			f_new[id0] = f0_new + eps
			f_new[id1] = f1_new - eps
			
			# Append to the output
			assert np.isclose(sum(f_new), 1.), "Molar fractions do not sum to unity : {}. Sum : {:.4f}".format(f_new, sum(f_new))
			fs.append(f_new)
			
	return np.array(fs)

def maximize_molar_fraction(func, n_elems, grid_step=0.05, step_threshold=0.005, func_args=()):
	'''
	Optimize ´func´ by first searching a regular grid with step size ´grid_step´.
	The maximum from the grid search is then used and optimized further on a grid,
	until the step size	has reached ´step_threshold´.
	Arguments to to ´func´ is passed via ´func_args´.
	'''
	
	# Get grid of molar fractions
	fs_grid = get_molar_fractions(grid_step, n_elems)
	
	# Convert to cartesian coordinates
	rs_grid = molar_fractions_to_cartesians(fs_grid)
	
	# Get function values on the grid
	ys_grid = func(rs_grid, *func_args)
	
	# Pick the point with the largest value on the grid
	idx_max = np.argmax(ys_grid)
	r_max = rs_grid[idx_max]
	f_max = cartesians_to_molar_fractions(r_max.reshape(1, -1))
	
	# Half the molar fraction step size
	step_size = grid_step / 2.
	
	# Optimize around the found maximum in ever decreasing steps
	while step_size > step_threshold:
		
		# Get molar fractions around the preliminary optimum
		fs_around = get_molar_fractions_around(f_max[0], step_size)
		
		# Convert to cartesian coordinates
		rs_around = molar_fractions_to_cartesians(fs_around)
		
		# Get function values of these coordinates
		ys_around = func(rs_around, *func_args)
		
		# Get index of the maximum value
		idx_max = np.argmax(ys_around)
		
		# Get cartesian coordinates of maximum
		r_max = rs_around[idx_max]
		
		# Convert maximum to molar fractions
		f_max = cartesians_to_molar_fractions(r_max.reshape(1, -1))
		
		# Half the molar fraction step size
		step_size /= 2.

	# Return molar fractions that maximize the function
	return f_max
	

# def maximize_molar_fraction_hypervolume(func, n_elems, grid_step=0.05, step_threshold=0.005, use_cartesian=True ,func_args=()):
# 	'''
# 	Optimize ´func´ by first searching a regular grid with step size ´grid_step´.
# 	The maximum from the grid search is then used and optimized further on a grid,
# 	until the step size	has reached ´step_threshold´.
# 	Arguments to to ´func´ is passed via ´func_args´.
# 	'''
	
# 	# Get grid of molar fractions
# 	fs_grid = get_molar_fractions(grid_step, n_elems)
	
# 	# Convert to cartesian coordinates
# 	if use_cartesian:
# 		rs_grid = molar_fractions_to_cartesians(fs_grid)
	
# 	# Get function values on the grid
# 	if use_cartesian:
# 		ys_grid = func(torch.tensor(rs_grid.reshape(rs_grid.shape[0],1,rs_grid.shape[1]),dtype=torch.double), *func_args).detach().numpy()
# 	else:
# 		ys_grid = func(torch.tensor(fs_grid.reshape(fs_grid.shape[0],1,fs_grid.shape[1]),dtype=torch.double), *func_args).detach().numpy()


# 	# Pick the point with the largest value on the grid
# 	y_max = np.max(ys_grid)
# 	idx_max = np.argmax(ys_grid)

# 	if use_cartesian:
# 		r_max = rs_grid[idx_max]
# 		f_max = cartesians_to_molar_fractions(r_max.reshape(1, -1))
# 	else:
# 		f_max = fs_grid[idx_max]
	
# 	# Half the molar fraction step size
# 	step_size = grid_step / 2.
	
# 	# Optimize around the found maximum in ever decreasing steps
# 	while step_size > step_threshold:
		
# 		# Get molar fractions around the preliminary optimum
# 		fs_around = get_molar_fractions_around(f_max.flatten(), step_size)
		
# 		# Convert to cartesian coordinates
# 		if use_cartesian:
# 			rs_around = molar_fractions_to_cartesians(fs_around)
		
# 		# Get function values of these coordinates
# 		if use_cartesian:
# 			ys_around = func(torch.tensor(rs_around.reshape(rs_around.shape[0],1,rs_around.shape[1]),dtype=torch.double), *func_args).detach().numpy()
# 		else:
# 			ys_around = func(torch.tensor(fs_around.reshape(fs_around.shape[0],1,fs_around.shape[1]),dtype=torch.double), *func_args).detach().numpy()

		
# 		# Check if molar fractions around are better than the preliminary optimimum
# 		y_max_around = np.max(ys_around)
# 		if y_max_around>y_max:
# 			y_max = y_max_around
# 			# Get index of the maximum value
# 			idx_max = np.argmax(ys_around)

# 			if use_cartesian:
# 				# Get cartesian coordinates of maximum
# 				r_max = rs_around[idx_max]
				
# 				# Convert maximum to molar fractions
# 				f_max = cartesians_to_molar_fractions(r_max.reshape(1, -1))
# 			else:
# 				f_max = fs_around[idx_max]

# 		# Half the molar fraction step size
# 		step_size /= 2.

# 	# Return molar fractions that maximize the function
# 	return f_max.flatten()

def make_triangle_ticks(ax, start, stop, tick, n, offset=(0., 0.),
						fontsize=12, ha='center', tick_labels=True):
	r = np.linspace(0, 1, n+1)
	x = start[0] * (1 - r) + stop[0] * r
	x = np.vstack((x, x + tick[0]))
	y = start[1] * (1 - r) + stop[1] * r
	y = np.vstack((y, y + tick[1]))
	ax.plot(x, y, 'black', lw=1., zorder=1)
	
	if tick_labels:
	
		# Add tick labels
		for xx, yy, rr in zip(x[0], y[0], r):
			ax.text(xx+offset[0], yy+offset[1], f'{rr*100.:.0f}',
					fontsize=fontsize, ha=ha)

def mscatter(x, y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    ax = ax or plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
	new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
		'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
		cmap(np.linspace(minval, maxval, n)))
	return new_cmap

def prepare_triangle_plot(ax, elems,edges=True,ticks=True):
	
	# Set the number of ticks to make
	n_ticks = 5
	tick_labels = True

	# Specify vertices as molar fractions
	fs_vertices = [[1., 0., 0.],
				   [0., 1., 0.],
				   [0., 0., 1.]]
	
	# Get cartesian coordinates of vertices
	xs_vertices, ys_vertices = molar_fractions_to_cartesians(fs_vertices).T
	
	# Get height of triangle
	h = 3**0.5/2
	
	# Define padding to put the vertex text neatly
	pad = [[-0.06, -0.06],
		   [ 0.06, -0.06],
		   [ 0.00,  0.08]]
	has = ['right', 'left', 'center']
	vas = ['top', 'top', 'bottom']

	# Make ticks and tick labels on the triangle axes
	left, right, top = np.concatenate((xs_vertices.reshape(-1,1), ys_vertices.reshape(-1,1)), axis=1)

	tick_size = 0.035
	bottom_ticks = 0.8264*tick_size * (right - top)
	right_ticks = 0.8264*tick_size * (top - left)
	left_ticks = 0.8264*tick_size * (left - right)

	# Set axis limits
	ax.set_xlim(-0.05, 1.05)
	ax.set_ylim(-0.05, h+0.05)
	if edges:
		# Plot triangle edges
		ax.plot([0., 0.5], [0., h], '-', color='black', zorder=1)
		ax.plot([0.5, 1.], [h, 0.], '-', color='black', zorder=1)
		ax.plot([0., 1.], [0., 0.], '-', color='black', zorder=1)
	
	# Remove spines
	for direction in ['right', 'left', 'top', 'bottom']:
		ax.spines[direction].set_visible(False)
	
	# Remove tick and tick labels
	ax.tick_params(which='both', bottom=False, labelbottom=False, left=False, labelleft=False)
	ax.set_aspect('equal')

	if ticks:		
		make_triangle_ticks(ax, right, left, bottom_ticks, n_ticks, offset=(0.03, -0.08), ha='center', tick_labels=tick_labels)
		make_triangle_ticks(ax, left, top, left_ticks, n_ticks, offset=(-0.03, -0.015), ha='right', tick_labels=tick_labels)
		make_triangle_ticks(ax, top, right, right_ticks, n_ticks, offset=(0.015, 0.02), ha='left', tick_labels=tick_labels)

		# Show axis labels (i.e. atomic percentages)
		ax.text(0.5, -0.14, f'{elems[0]} content (at.%)', rotation=0., fontsize=12, ha='center', va='center')
		ax.text(0.9, 0.5, f'{elems[1]} content (at.%)', rotation=-60., fontsize=12, ha='center', va='center')
		ax.text(0.1, 0.5, f'{elems[2]} content (at.%)', rotation=60., fontsize=12, ha='center', va='center')
	else:
		pad = np.zeros((3,2))
	# Show the chemical symbol as text at each vertex
	for idx, (x, y, (dx, dy)) in enumerate(zip(xs_vertices, ys_vertices, pad)):
		ax.text(x+dx, y+dy, s=elems[idx], fontsize=14, ha=has[idx], va=vas[idx])
	
	return ax

			  
def prepare_tetrahedron_plot(ax, elems):
	
	# Set the number of ticks to make
	n_ticks = 5
	tick_labels = True

	# Specify vertices as molar fractions
	fs_vertices = [[1., 0., 0., 0.],
				   [0., 1., 0., 0.],
				   [0., 0., 1., 0.],
				   [0., 0., 0., 1.]]
	
	# Get cartesian coordinates of vertices
	rs_vertices = molar_fractions_to_cartesians(fs_vertices)

	# Set axis limits
	ax.set_xlim(-0.05, np.max(rs_vertices[:, 0]) + 0.05)
	ax.set_ylim(-0.05, np.max(rs_vertices[:, 1]) + 0.05)
	ax.set_zlim(-0.05, np.max(rs_vertices[:, 2]) + 0.05)

	# Get all combinations of pairs of vertices
	for edge in it.combinations(rs_vertices, 2):
	
		# Plot simplex edge
		ax.plot(*np.array(edge).T, ls='solid', color='black')
	
	# Remove spines
	for direction in ['right', 'left', 'top', 'bottom']:
		ax.spines[direction].set_visible(False)
	
	# Remove spines etc.
	ax.set_axis_off()
	
	# Define padding to put the vertex text neatly
	pad = [[-0.03, -0.03, 0.00],
		   [ 0.03, -0.03, 0.00],
		   [ 0.03,  0.03, 0.00],
		   [ 0.00,  0.00, 0.06]]
	has = ['right', 'left', 'center', 'center']
	vas = ['top', 'top', 'bottom', 'center']
	
	# Show the chemical symbol as text at each vertex
	for idx, (r, dr) in enumerate(zip(rs_vertices, pad)):
		ax.text(*(r+dr).T, s=elems[idx], fontsize=14, ha=has[idx], va=vas[idx])
	
	return ax




def make_ternary_plot(grid,target_func,elements,ax=None,contour_levels=15,vmin=0,vmax=None,colorbar=False,colormap='viridis',minval=0.2,maxval=1.0,cbar_ticks=None,cbar_label=None):
	# Define color map of plot
	cmap = truncate_colormap(plt.get_cmap(colormap), minval=minval, maxval=maxval, n=100)

	if ax is None:
		# Make figure for plotting Gaussian process mean, uncertainty, and acquisition function
		fig, ax = plt.subplots(figsize=(4,4),dpi=400)
		return_fig=True
	else:
		return_fig=False

	# Apply shared settings to pseudo-ternary plot
	prepare_triangle_plot(ax, elems=elements)

	# Plot surrogate/uncertainty/acquisition function as a contour plot
	plot_kwargs = dict(cmap=cmap, levels=contour_levels, zorder=0,vmin=vmin,vmax=vmax)

	contour=ax.tricontourf(*grid, target_func, **plot_kwargs)

	if colorbar:
		plt.colorbar(contour,ax=ax,shrink=0.5,anchor=(0.0,0.85),ticks=cbar_ticks,label=cbar_label)

	if return_fig:
		return fig, ax
	else:
		return ax