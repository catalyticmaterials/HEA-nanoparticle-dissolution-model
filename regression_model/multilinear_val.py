from utilities import metals,U_standard
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
import matplotlib.pyplot as plt
from parity_plot import partity_plot

n_metals = len(metals)
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
metal_U = metal_data[:,n_metals:].flatten()

features = np.vstack((features,features_metal))
dE = np.append(dE,metal_dE)

metal_feature = np.vstack((metal_feature,np.vstack([np.eye(n_metals)]*7)))



# regressor = MultiLinearRegressor()
regressor = LinearRegression(fit_intercept=False)
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
from scipy.optimize import minimize

# constraints = LinearConstraint(A,ub=np.zeros(6),keep_feasible=True)

loo = LeaveOneOut()
preds = []
for train_index,test_index in loo.split(features):
    # print(test_index)
    X_train = features[train_index]
    X_test = features[test_index]

    # y_train = U[train_index]
    y_train = dE[train_index]

    metal_idx = np.where(metal_feature[test_index][0]==1)[0][0]


    mask = metal_feature[train_index,metal_idx]==1


    def metal_constrain(params):
        return params[7+metal_idx]

    # Define the constraint dictionary in the form required by 'minimize'
    constraints = [{'type': 'ineq', 'fun': combined_constraint},
                   {'type': 'eq','fun':eq_constrain},
                   {'type': 'eq','fun':metal_constrain}]
    result = minimize(cost_function, initial_params.copy(), args=(X_train[mask], y_train[mask]),constraints=constraints,tol=1e-12,method='SLSQP')
    
    params = result.x


    
    preds.append(np.dot(X_test,params)[0])


# mask = dE[:100]>-1.5
# print(np.mean(np.abs(np.array(preds)-dE)))


partity_plot(dE,preds,metal_feature,unit='eV')

# plt.tight_layout()
plt.savefig('loocv_parity_eV.png',dpi=600,bbox_inches='tight')


preds = np.array(preds)

preds_U = np.zeros(len(preds))
U =np.zeros_like(preds_U)
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds_U[mask] = np.min(preds[mask]/U_standard[metal]['n'] + U_standard[metal]['U'],axis=0)
    # preds_U[mask] = U_diss(preds[mask],metal,1e-6)
    U[mask] = np.min(dE[mask]/U_standard[metal]['n'] + U_standard[metal]['U'],axis=0)
    # U[mask] = U_diss(dE[mask],metal,1e-6)

partity_plot(U,preds_U,metal_feature,unit='V')

# plt.tight_layout()
plt.savefig('loocv_parity_V.png',dpi=600,bbox_inches='tight')


