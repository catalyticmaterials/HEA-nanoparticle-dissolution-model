import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF,ConstantKernel,WhiteKernel
from utilities import metals

gpr = GPR(kernel=ConstantKernel(1.0) * RBF(1.0) + WhiteKernel(noise_level=0.05**2,noise_level_bounds=(1e-9,0.1)),n_restarts_optimizer=25,alpha=0)


train_data=np.loadtxt('full_space_sampling.csv',skiprows=1,delimiter=',')
X_train = train_data[:,:-1]
y_train = train_data[:,-1]

gpr.fit(X_train,y_train)

for quinary_metals,name in zip((['Ag','Au','Cu','Pd','Pt'],['Ir','Pd','Pt','Rh','Ru']),('GPGM','PGM')):

    metal_mask = np.array([metal in quinary_metals for metal in metals])


    test_data = np.loadtxt(f'{name}_test_sampling.csv',skiprows=1,delimiter=',')

    X_test = np.zeros((100,8))
    X_test[:,metal_mask] = test_data[:,:-1]
    y_test = test_data[:,-1]

    y_pred = gpr.predict(X_test)

    E = y_pred - y_test

    ME = np.mean(E)
    MAE = np.mean(np.abs(E))
    RMSE = np.sqrt(np.mean(E**2))

    print(name)
    print('Full space GPR')
    print(f'ME: {ME:.3f}\nMAE: {MAE:.3f}\nRMSE: {RMSE:.3f}')



    train_data = np.loadtxt(f'{name}_train_sampling.csv',skiprows=1,delimiter=',')
    gpr_sub = GPR(kernel=ConstantKernel(1.0) * RBF(1.0) + WhiteKernel(noise_level=0.05**2,noise_level_bounds=(1e-9,0.1)),n_restarts_optimizer=25,alpha=0)
    X_train = train_data[:,:-1]
    y_train = train_data[:,-1]
    gpr_sub.fit(X_train,y_train)

    y_pred = gpr_sub.predict(X_test[:,metal_mask])

    E = y_pred - y_test

    ME = np.mean(E)
    MAE = np.mean(np.abs(E))
    RMSE = np.sqrt(np.mean(E**2))

    print('Subspace GPR')
    print(f'ME: {ME:.3f}\nMAE: {MAE:.3f}\nRMSE: {RMSE:.3f}\n')

