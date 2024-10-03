import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

data = np.loadtxt('Sd_variance.csv',skiprows=1,delimiter=',')
mf = data[:,:-2]
Sd_mean = data[:,-2]
Sd_std = data[:,-1]

# idx = np.argsort(Sd_std)[-4:]
# print(mf[idx])
# print(Sd_mean[idx])
# print(Sd_std[idx])

# idx = np.argsort(Sd_std)[:4]
# print(mf[idx])
# print(Sd_mean[idx])
# print(Sd_std[idx])

# plt.hist(Sd_std,bins=15)
# plt.show()


from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler

for y in (Sd_mean,Sd_std):
    # Correlation
    # Calculate correlations
    correlations = []
    p_values = []

    for i in range(mf.shape[1]):
        corr, p_value = pearsonr(mf[:, i], y)
        correlations.append(corr)
        p_values.append(p_value)

    # Output the results
    for i, (corr, p_value) in enumerate(zip(correlations, p_values)):
        print(f"Feature {i+1}: Correlation = {corr:.3f}, p-value = {p_value:.3f}")






    # Standardize the features
    scaler = StandardScaler()
    X_standardized = scaler.fit_transform(mf)

    # Apply PCA
    pca = PCA(2)
    X_pca = pca.fit_transform(X_standardized)



    correlations = []
    p_values = []

    for i in range(X_pca.shape[1]):
        corr, p_value = pearsonr(X_pca[:, i], y)
        correlations.append(corr)
        p_values.append(p_value)


    # Output the correlation results
    for i, (corr, p_value) in enumerate(zip(correlations, p_values)):
        print(f"Principal Component {i+1}: Correlation = {corr:.3f}, p-value = {p_value:.3f}")

    # Output the explained variance of each principal component
    explained_variance = pca.explained_variance_ratio_
    for i, var in enumerate(explained_variance):
        print(f"Principal Component {i+1}: Explained Variance = {var:.3f}")


    print(pca.components_)
    print()







# Correlation
    # Calculate correlations
    correlations = []
    p_values = []

    for i in range(mf.shape[1]):
        corr, p_value = pearsonr(mf[:, i], y)
        correlations.append(corr)
        p_values.append(p_value)

    # Output the results
    for i, (corr, p_value) in enumerate(zip(correlations, p_values)):
        print(f"Feature {i+1}: Correlation = {corr:.3f}, p-value = {p_value:.3f}")





train_data=np.loadtxt('full_space_sampling.csv',skiprows=1,delimiter=',')
mf = train_data[:,:-1]
y = train_data[:,-1]


# Standardize the features
scaler = StandardScaler()
X_standardized = scaler.fit_transform(mf)

# Apply PCA
pca = PCA(2)
X_pca = pca.fit_transform(X_standardized)



correlations = []
p_values = []

for i in range(X_pca.shape[1]):
    corr, p_value = pearsonr(X_pca[:, i], y)
    correlations.append(corr)
    p_values.append(p_value)


# Output the correlation results
for i, (corr, p_value) in enumerate(zip(correlations, p_values)):
    print(f"Principal Component {i+1}: Correlation = {corr:.3f}, p-value = {p_value:.3f}")

# Output the explained variance of each principal component
explained_variance = pca.explained_variance_ratio_
for i, var in enumerate(explained_variance):
    print(f"Principal Component {i+1}: Explained Variance = {var:.3f}")


print(pca.components_)
print()



