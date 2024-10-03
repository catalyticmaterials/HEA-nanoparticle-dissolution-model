from.colors import metal_colors
import pickle
import numpy as np

# metals = ['Ir','Pd','Pt','Rh','Ru']
metals = ['Ag','Au','Cu','Ir','Pd','Pt','Rh','Ru']

# Standard reduction potentials from CRC handbook 97th edition (2016)
# U_standard = { 
#     'Ag': {'n':1, 'U':0.80},
#     'Au': {'n':3,'U':1.50},
#     'Cu': {'n':2,'U':0.34},
#     'Ir':{'n':3,'U':1.16},
#     'Pd':{'n':2,'U':0.95},
#     'Pt':{'n':2,'U':1.18}, 
#     'Rh':{'n':1,'U':0.60},
#     'Ru':{'n':2,'U':0.46}
# }

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


# Standard reduction potentials with Cl from CRC handbook 97th edition (2016)
U_standard_Cl = {
    'Ag': {'n':1, 'U':0.22233},
    'Au': {'n':3,'U':1.002},
    'Ir':{'n':3,'U':0.77},
    'Pd':{'n':2,'U':0.591},
    'Pt':{'n':2,'U':0.755}, 
    'Rh':{'n':3,'U':0.431},
    'Ru':{'n':2,'U':0.455},
}

def load_regressor(model_path):
    
    with open(model_path,'rb') as input:
        regressor = pickle.load(input)

    return regressor


kB=8.617333262e-5
T=298.15

def U_diss(dE,metal,c):
    return np.min(dE/U_standard[metal]['n'] + U_standard[metal]['U'] + kB*T/U_standard[metal]['n'] * np.log(c),axis=0)