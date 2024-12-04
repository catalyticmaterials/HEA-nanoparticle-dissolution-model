from.colors import metal_colors
import pickle
import numpy as np

metals = ['Ag','Au','Cu','Ir','Pd','Pt','Rh','Ru']

# Standard reduction potentials from CRC handbook 97th edition (2016)
U_standard = { 
    'Ag': {'n':np.array([[1]]), 'U':np.array([[0.80]])},
    'Au': {'n':np.array([[1],[3]]),'U':np.array([[1.69],[1.50]])},
    'Cu': {'n':np.array([[1],[2]]),'U':np.array([[0.52],[0.34]])},
    'Ir': {'n':np.array([[3]]),'U':np.array([[1.16]])},
    'Pd': {'n':np.array([[2]]),'U':np.array([[0.95]])},
    'Pt': {'n':np.array([[2]]),'U':np.array([[1.18]])}, 
    'Rh': {'n':np.array([[1],[3]]),'U':np.array([[0.60],[0.76]])},
    'Ru': {'n':np.array([[2]]),'U':np.array([[0.46]])}
}



def load_regressor(model_path):
    
    with open(model_path,'rb') as input:
        regressor = pickle.load(input)

    return regressor


kB=8.617333262e-5
T=298.15

def U_diss(dE,metal,c):
    return np.min(dE/U_standard[metal]['n'] + U_standard[metal]['U'] + kB*T/U_standard[metal]['n'] * np.log(c),axis=0)