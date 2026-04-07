import numpy as np

def makeNoiseMatrix(n, randomVector, constant):
	I=np.zeros((n,n))
	for i in range(0,n):
		for j in range(0,n):
			if i==j:
				I[i][j] = randomVector[i]*constant
			else:
				I[i][j] = 0
	return I

def subVector(x1, x2, size):
	
	outVector=np.zeros(size)
	for i in range (0, size):
		outVector[i]=x1[i]-x2[i]
	
	return outVector

def addMatrix(x1, x2, sizeY, sizeX): 
	
	outMatrix=np.zeros((sizeY,sizeX))
	for	i in range(0,sizeY):
		for j in range (0,sizeX):
			outMatrix[i][j]= x1[i][j]+x2[i][j]
	return outMatrix

def logDet(a, n):
	#log(det(A)) = 2*sum(log(vecdiag(G))) 
	#https://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix.html
	chol = np.linalg.cholesky(a)
	result = 0
	
	for i in range(0,n):
		result += np.log(chol[i][i])

	result = result * 2
	return result

def traceMatrix(Matrix, size):
	result = 0
	for i in range(0,size):
		result += Matrix[i][i]
	return result

def multVector(Matrix1, Vector, dim_1, dim_2) :
	vectorDouble=np.zeros(dim_1)
 
	for i in range (0,dim_1):
		suma = 0
		for j in range (0, dim_2):
			suma += Matrix1[i][ j] * Vector[j]
		vectorDouble[i] = suma

	return vectorDouble
