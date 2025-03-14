from utilities.linear_model import MultiLinearRegressor
from utilities import metals, U_diss
from parity_plot import partity_plot
import matplotlib.pyplot as plt
import numpy as np
import pickle

n_metals = len(metals)

regressor = MultiLinearRegressor()

data = np.loadtxt('../DFT_calculations/hea_data.csv',delimiter=',',skiprows=1)

features_ = data[:,:-3]

dE = data[:,-3]
cn = np.asarray(np.sum(features_[:,n_metals:],axis=1),dtype=int)


cn_feature = np.zeros((len(cn),7))
cn_feature[cn==3,0] = 1
cn_feature[cn==4,1] = 1
cn_feature[cn==5,2] = 1
cn_feature[cn==6,3] = 1
cn_feature[cn==7,4] = 1
cn_feature[cn==8,5] = 1
cn_feature[cn==9,6] = 1

metal_feature = features_[:,:n_metals]
features = np.hstack((cn_feature,features_[:,n_metals:]/cn.reshape(-1,1)))

# features for metal data
I_cn = np.eye(7)
cn3 = np.hstack(([I_cn[0]]*n_metals,np.eye(n_metals)))
cn4 = np.hstack(([I_cn[1]]*n_metals,np.eye(n_metals)))
cn5 = np.hstack(([I_cn[2]]*n_metals,np.eye(n_metals)))
cn6 = np.hstack(([I_cn[3]]*n_metals,np.eye(n_metals)))
cn7 = np.hstack(([I_cn[4]]*n_metals,np.eye(n_metals)))
cn8 = np.hstack(([I_cn[5]]*n_metals,np.eye(n_metals)))
cn9 = np.hstack(([I_cn[6]]*n_metals,np.eye(n_metals)))

features_metal = np.vstack((cn3,cn4,cn5,cn6,cn7,cn8,cn9))

metal_data = np.loadtxt('../DFT_calculations/metals_data.csv',delimiter=',')
metal_dE = metal_data[:,:n_metals].flatten()


features = np.vstack((features,features_metal))
dE = np.append(dE,metal_dE)
metal_feature = np.vstack((metal_feature,np.vstack([np.eye(n_metals)]*7)))


# Cost function: sum of squared residuals
def cost_function(params, X, y):
    predictions = X @ params
    return np.sum((y - predictions) ** 2)

# Combined constraint function: ensures that params[i+1] >= params[i]
def combined_constraint(params):
    return np.diff(params[:7])  # Returns an array of differences: params[i+1] - params[i]

def eq_constrain(params):
    return params[3]



# Initial guess for the parameters
initial_params = np.append(np.arange(-3,4)/3,np.zeros(8))

for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    def metal_constrain(params):
        return params[7+i]

    # Define the constraint dictionary in the form required by 'minimize'
    constraints = [{'type': 'ineq', 'fun': combined_constraint},
                   {'type': 'eq','fun':eq_constrain},
                   {'type': 'eq','fun':metal_constrain}]
    result = regressor.fit(metal,features[mask],dE[mask],cost_function,constraints,initial_params)
    
    


# with open('../utilities/AgAuCuIrPdPtRhRu_multilinear.regressor','wb') as out:
# 	pickle.dump(regressor,out)
      


preds = np.zeros(len(dE))
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds[mask] = regressor.predict(metal,features[mask])


partity_plot(dE,preds,metal_feature)

plt.savefig('parity_plots/train_parity_eV.png',dpi=600,bbox_inches='tight')


preds = np.array(preds)

preds_U = np.zeros(len(preds))
U = np.zeros_like(preds_U)
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds_U[mask] = U_diss(preds[mask],metal,1e-6)
    U[mask] = U_diss(dE[mask],metal,1e-6)

partity_plot(U,preds_U,metal_feature,unit='V')

plt.savefig('parity_plots/train_parity_V.png',dpi=600,bbox_inches='tight')


