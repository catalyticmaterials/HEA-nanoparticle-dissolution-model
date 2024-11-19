from scipy.optimize import minimize

class MultiLinearRegressor():
	'Class for doing multilinear regression'
	
	def __init__(self):
		self.parameters = {}
		
	def fit(self, metal, X, y,loss_func,constraints,x0):
		
		'Train regressor on X and y for the specified metal'
		self.parameters[metal] = minimize(loss_func, x0, args=(X, y),constraints=constraints,tol=1e-12,method='SLSQP').x
		
	def predict(self, metal, X):
		'Predict on X for the specified ensemble'
		return self.parameters[metal] @ X.T
	
	