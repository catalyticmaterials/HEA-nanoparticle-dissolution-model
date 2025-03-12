from utilities import metals,U_diss
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import confusion_matrix, precision_score, recall_score, accuracy_score
import matplotlib.pyplot as plt
from parity_plot import partity_plot
from scipy.optimize import minimize


n_metals = len(metals)
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
metal_U = metal_data[:,n_metals:].flatten()

features = np.vstack((features,features_metal))
dE = np.append(dE,metal_dE)

metal_feature = np.vstack((metal_feature,np.vstack([np.eye(n_metals)]*7)))




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




loo = LeaveOneOut()
preds = []
for train_index,test_index in loo.split(features):

    X_train = features[train_index]
    X_test = features[test_index]

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


partity_plot(dE,preds,metal_feature,unit='eV')


plt.savefig('parity_plots/loocv_parity_eV.png',dpi=600,bbox_inches='tight')


preds = np.array(preds)

preds_U = np.zeros(len(preds))
U =np.zeros_like(preds_U)
for i,metal in enumerate(metals):
    mask = metal_feature[:,i]==1
    preds_U[mask] = U_diss(preds[mask],metal,1e-6)
    U[mask] = U_diss(dE[mask],metal,1e-6)

partity_plot(U,preds_U,metal_feature,unit='V')

plt.savefig('parity_plots/loocv_parity_V.png',dpi=600,bbox_inches='tight')

dE_mask = dE>-0.5
print('MAE for dE above -0.5 eV:',np.mean(np.abs(dE[dE_mask]-preds[dE_mask])))
for i,metal in enumerate(metals):
    metal_mask = metal_feature[:,i]==1
    mask = metal_mask*dE_mask
    print(metal,np.mean(np.abs(dE[mask]-preds[mask])))


U_mask = ((U>=0.7)*(U<=0.9))*((preds_U>=0.7)*(preds_U<=0.9))
print('MAE for calc and predicted Udiss between 0.7 and 0.9 V:',np.mean(np.abs(U[U_mask]-preds_U[U_mask])))
for i,metal in enumerate(metals):
    metal_mask = metal_feature[:,i]==1
    mask = metal_mask*U_mask
    print(metal,np.mean(np.abs(U[mask]-preds_U[mask])))




# Plot confusion matrix
U_target = 0.8
dissolve_true = U<U_target
dissolve_pred = preds_U<U_target




# Compute confusion matrix
cm = confusion_matrix(dissolve_true, dissolve_pred)
tn, fp, fn, tp = cm.ravel()
cm = np.array([[tp,fn],[fp,tn]])
# Compute precision, recall and accuracy
precision = precision_score(dissolve_true, dissolve_pred)
recall = recall_score(dissolve_true, dissolve_pred)
accuracy = accuracy_score(dissolve_true, dissolve_pred)

# Plot confusion matrix using matplotlib cmap
plt.figure(figsize=(4,4))
plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
plt.title(f'$U_{{diss}}$ versus 0.8 V\nPrecision: {precision:.2f}, Recall: {recall:.2f}, Accuracy: {accuracy:.2f}',fontsize='medium')
tick_marks = (0,1)
tick_labels = ['Dissolved','Stable']
plt.xticks(tick_marks, tick_labels)
plt.yticks(tick_marks, tick_labels,rotation=90,va='center')

# Annotate the confusion matrix
thresh = cm.max() / 2.
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        plt.text(j, i, format(cm[i, j], 'd'),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

plt.ylabel('DFT')
plt.xlabel('Predicted')
plt.tight_layout()
plt.savefig('parity_plots/Confusion_matrix.png',dpi=600)