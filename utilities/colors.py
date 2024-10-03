from matplotlib.colors import to_rgb
import numpy as np

# Specify colors of metals
metal_colors = dict(Ag='silver',
					Au='gold',
					Cu='darkorange',
                    Ir='darkslateblue',
					Pd='royalblue',
					Pt='green',
                    Rh='teal',
					Ru='darkred')

def hollow_site_color(site_metals):

    return np.mean(np.array([to_rgb(metal_colors[metal]) for metal in site_metals]),axis=0)

def alloy_color(metals,compositions):
    return np.sum(np.array([np.array(to_rgb(metal_colors[metal]))*x for metal,x in zip(metals,compositions)]),axis=0)