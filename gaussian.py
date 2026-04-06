import numpy as np
from mathFunctions import *
from kernelFunctions import *

def KernelMV(xKernel, x2Kernel, size_x, size_x2, Dimension, k_param, type):
	cov=np.zeros((size_x,size_x2))
	for i in range(0,size_x):
		for j in range(0,size_x2):
			if type==0:
				cov[i][j] = maternKernel52(xKernel, x2Kernel, i, j, Dimension, k_param)
			if type==1:
				cov[i][j] = exponentialKernel(xKernel, x2Kernel, i, j, Dimension, k_param)
			if type==2:
				cov[i][j] = maternKernel32(xKernel, x2Kernel, i, j, Dimension, k_param)
	return cov

def marginalLikelihood(params,n_points, xtrain, Dimension,  ytrain, type, noiseLevel):
	randomVec =np.random.rand(n_points)
	noise = makeNoiseMatrix(n_points, randomVec, noiseLevel)
	#print(noise)
	K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type)
	#print(K)
	A = addMatrix(K, noise, n_points, n_points)
	#print(A)
	Ainv = np.linalg.inv(A) 
	Ainvy=np.matmul(Ainv, ytrain)
	yt_Ainvy=np.dot(ytrain,Ainvy)
	
	logDeterminant = logDet(A, n_points)
	marginal_likelihood = (-0.5*yt_Ainvy - 0.5*(logDeterminant) - 0.5 *np.log(2 * np.pi)*n_points)

	return marginal_likelihood
	
def marginalLikelihoodNeg(params,n_points, xtrain, Dimension,  ytrain, type, noiseLevel):
	randomVec =np.random.rand(n_points)
	noise = makeNoiseMatrix(n_points, randomVec, noiseLevel)
	#print(noise)
	K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type)
	#print(K)
	A = addMatrix(K, noise, n_points, n_points)
	#print(A)
	Ainv = np.linalg.inv(A) 
	Ainvy=np.matmul(Ainv, ytrain)
	yt_Ainvy=np.dot(ytrain,Ainvy)
	
	logDeterminant = logDet(A, n_points)
	marginal_likelihood = (-0.5*yt_Ainvy - 0.5*(logDeterminant) - 0.5 *np.log(2 * np.pi)*n_points)

	return -marginal_likelihood

def diffKernelMV(xKernel, x2Kernel, size_x, size_x2, Dimension, k_param, type, kernelIndex) : 
	cov=np.zeros((size_x,size_x2))
	for i in range(0,size_x):
		for j in range(0,size_x2):
			if type==0:
				cov[i][j] = diffMaternKernel52(xKernel, x2Kernel, i, j, Dimension, k_param, kernelIndex)
			if type==1:
				cov[i][j] = diffExpKernel(xKernel, x2Kernel, i, j, Dimension, k_param, kernelIndex)
			if type==2:
				cov[i][j] = diffMaternKernel32(xKernel, x2Kernel, i, j, Dimension, k_param, kernelIndex)
	return cov

def marginalLikelihoodDiff(params, n_points, xtrain, Dimension,  ytrain, type, noiseLevel):
	randomVec =np.random.rand(n_points)
	noise = makeNoiseMatrix(n_points, randomVec, noiseLevel)
	#print(noise)
	K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type)
	#print(K)
	A = addMatrix(K, noise, n_points, n_points)
	#print(A)
	Ainv = np.linalg.inv(A) 
	Ainvy=np.matmul(Ainv, ytrain)
	yt_Ainvy=np.dot(ytrain,Ainvy)
	
	logDeterminant = logDet(A, n_points)
	marginal_likelihood = (-0.5*yt_Ainvy - 0.5*(logDeterminant) - 0.5 *np.log(2 * np.pi)*n_points)
	f = -(marginal_likelihood)
	
	
	g=np.zeros(Dimension+1)
	
	for indexTheta in range(0, Dimension+1):
		pK_Theta = diffKernelMV(xtrain, xtrain,n_points, n_points, Dimension, params, type, indexTheta)
		tempTraceTheta = np.matmul(Ainv, pK_Theta)
		traceTheta = traceMatrix(tempTraceTheta, n_points)
		pTheta = multVector(pK_Theta, Ainvy, n_points, n_points)
		AinvPTheta = multVector(Ainv, pTheta, n_points, n_points)
		yAinvPTheta = np.dot(ytrain, AinvPTheta)

		g[indexTheta] = -(0.5 * yAinvPTheta - 0.5 * traceTheta);
		
	return g

def getAinv(n_points, xtrain, Dimension, params, ytrain, type, noiseLevel):
	randomVec =np.random.rand(n_points)
	noise = makeNoiseMatrix(n_points, randomVec, noiseLevel)
	#print(noise)
	K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type)
	
	A = addMatrix(K, noise, n_points, n_points)
	#print(A)
	Ainv = np.linalg.inv(A) 

	return Ainv
	
def acqFunctionUCB(testpoint, Ainv, wi, n_x, dimension, params, xtrain, kappa, type, bestSolution):

	#print("testpoint %s"%testpoint)
	testpointMat=np.zeros((1, dimension))

	for i in range(0,dimension):
		testpointMat[0][i]=testpoint[i]

	#Mean*******************************************************************************
	testpoint_kernel = KernelMV(xtrain, testpointMat, n_x, 1, dimension, params, type)

	transtestpoint_kernel = testpoint_kernel.transpose()

	mean_predTemp = multVector(transtestpoint_kernel, wi, 1, n_x)

	mean_pred = mean_predTemp[0]


	#Var*******************************************************************************
	Ainvk = np.matmul(Ainv, testpoint_kernel)

	ktAinvk = np.matmul(transtestpoint_kernel, Ainvk)

	kxx= KernelMV(testpointMat, testpointMat, 1, 1, dimension, params, type)
	
	var_pred = ((kxx[0][0] - ktAinvk[0][0]))

	#	/*printf("mean \t %0.8f \n", mean_pred);
	#	printf("var \t %0.8f \n", var_pred);*/

	result = (mean_pred + kappa*np.sqrt(var_pred));

	for i in range(0, dimension):
		if testpoint[i] > 1:
			result = 1e20
		if testpoint[i] < 0:
			result = 1e20

	return result
