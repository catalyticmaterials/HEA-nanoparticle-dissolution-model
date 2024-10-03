import numpy as np

cols = [0,1,2,4,5]

data3 = np.loadtxt('metals_data.csv',delimiter=',',skiprows=1,usecols=cols)

data4 = np.loadtxt('4x4/metals_data.csv',delimiter=',',skiprows=1,usecols=cols)


MAE = np.mean(np.abs(data3-data4))

print(MAE)
print(np.mean(data3-data4))

