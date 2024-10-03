import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from utilities.stability import metals
from ase.db import connect

dE = np.loadtxt('metals_data.csv',delimiter=',',skiprows=1,usecols=(0,1,2,3,4)).T
cn = np.arange(3,10)
# print(dE)
slopes = []
for i,metal in enumerate(metals):
    slopes.append(linregress(cn,dE[i]).slope)

with connect('bulks.db') as db:
    E_bulk = [row.energy for row in db.select()]

surfaces = ['T111','T100','edge','kink']
E_surfaces = {}
with connect('metals_dissolution_out.db') as db:
    for surface in surfaces:
        E_surface = []
        for i,row in enumerate(db.select(surface=surface,defect='none')):
            E_total = row.energy

            atoms = db.get_atoms(row.id)
            n=len(atoms)

            a,b = atoms.cell.cellpar()[:2]
            A = a*b

            E_surface.append((E_total - n*E_bulk[i])/(2*A))
        E_surfaces[surface] = E_surface



fig,ax = plt.subplots()
for surface in surfaces:
    ax.scatter(E_surfaces[surface],slopes,label=surface)

print(E_surfaces)

ax.legend()
plt.show()