from utilities.linear_model import MultiLinearRegressor
from utilities import metals,U_standard
from parity_plot import partity_plot
import matplotlib.pyplot as plt
import numpy as np
import pickle

n_metals = len(metals)

regressor = MultiLinearRegressor()

data = np.loadtxt('../DFT_calculations/hea_data.csv',delimiter=',',skiprows=1)

features_ = data[:,:-3]

# cn = data[:,-2]
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
features = np.zeros((len(features_),8*7))
for i,feature in enumerate(features_):
    features[i,8*(cn[i]-3):8*(cn[i]-3)+8] = feature[n_metals:]

# features = np.hstack((cn_feature,features_[:,n_metals:]/cn.reshape(-1,1)))
features = np.hstack((cn_feature,features))


# features for metal data
I_cn = np.eye(7)
I_metals = np.eye(56)

cn3 = np.hstack(([I_cn[0]]*n_metals,I_metals[0:8]*3))
cn4 = np.hstack(([I_cn[1]]*n_metals,I_metals[8:16]*4))
cn5 = np.hstack(([I_cn[2]]*n_metals,I_metals[16:24]*5))
cn6 = np.hstack(([I_cn[3]]*n_metals,I_metals[24:32]*6))
cn7 = np.hstack(([I_cn[4]]*n_metals,I_metals[32:40]*7))
cn8 = np.hstack(([I_cn[5]]*n_metals,I_metals[40:48]*8))
cn9 = np.hstack(([I_cn[6]]*n_metals,I_metals[48:56]*9))

features_metal = np.vstack((cn3,cn4,cn5,cn6,cn7,cn8,cn9))

metal_data = np.loadtxt('../DFT_calculations/metals_data.csv',delimiter=',')
metal_dE = metal_data[:,:n_metals].flatten()
metal_U = metal_data[:,n_metals:].flatten()

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
initial_params = np.append(np.arange(-3,4)/3,np.zeros(56))

for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    def metal_constrain(params):
        return [params[7+ j*8 + i] for j in range(7)]

    # Define the constraint dictionary in the form required by 'minimize'
    constraints = [{'type': 'ineq', 'fun': combined_constraint},
                   {'type': 'eq','fun':eq_constrain},
                   {'type': 'eq','fun':metal_constrain}]
    result = regressor.fit(metal,features[mask],dE[mask],cost_function,constraints,initial_params)
    
    

print(np.round(regressor.parameters['Pt'],decimals=3))
print(np.round(regressor.parameters['Au'],decimals=3))

Pt1 = regressor.predict('Pt',np.array([0,1,0,0,0,0,0] + [0]*8 + [0,2,0,0,0,1,0,1] + [0]*40))
Pt2 = regressor.predict('Pt',np.array([0,0,1,0,0,0,0] + [0]*16 +[0,3,0,0,0,1,0,1] + [0]*32))
Au1 = regressor.predict('Au',np.array([0,1,0,0,0,0,0] + [0]*8 + [0,2,0,0,1,0,0,1] + [0]*40))
Au2 = regressor.predict('Au',np.array([0,0,1,0,0,0,0] + [0]*16 +[0,2,0,0,1,1,0,1] + [0]*32))
# Au1 = regressor.predict('Au',np.array([0,0,0,0,1,0,0] + [0]*32 + [0,2,0,0,1,0,0,1] + [0]*16))
# Au2 = regressor.predict('Au',np.array([0,0,0,0,0,1,0] + [0]*40 +[0,2,0,0,1,1,0,1] + [0]*8))

print(Pt2-Pt1)
print(Au2-Au1)
stop
# with open('../utilities/AgAuCuIrPdPtRhRu_multilinear.regressor','wb') as out:
# 	pickle.dump(regressor,out)
      


preds = np.zeros(len(dE))
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds[mask] = regressor.predict(metal,features[mask])


partity_plot(dE,preds,metal_feature)

plt.tight_layout()
plt.savefig('train_parity_eV_copy2.png',dpi=600)


preds = np.array(preds)

preds_U = np.zeros(len(preds))
U = np.zeros_like(preds_U)
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds_U[mask] = np.min(preds[mask]/U_standard[metal]['n'] + U_standard[metal]['U'],axis=0)

    U[mask] = np.min(dE[mask]/U_standard[metal]['n'] + U_standard[metal]['U'],axis=0)


partity_plot(U,preds_U,metal_feature,unit='V')

plt.tight_layout()
plt.savefig('train_parity_V_copy2.png',dpi=600)


