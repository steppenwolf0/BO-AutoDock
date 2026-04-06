import numpy as np
from mathFunctions import *
import math 

def maternKernel52(xKernel, yKernel, indexI, indexJ, Dimension,  k_param):
	#This is the Matern 5/2 Kernel
	# https://nl.mathworks.com/help/stats/kernel-covariance-function-options.html
	xI = np.zeros(Dimension)
	xJ = np.zeros(Dimension)
	
	for i in range(0, Dimension):
		xI[i] = xKernel[indexI][i]
		xJ[i] = yKernel[indexJ][i]

	temp = subVector(xI, xJ, Dimension)

	r = 0
	for i in range(0, Dimension):
		x1 = temp[i]
		t1 = k_param[i]
		r += (x1 * x1) / (np.log(t1) * np.log(t1))

	r = np.sqrt(r)
	t0 =  k_param[Dimension];

	factor1 = 1 + np.sqrt(5) * r + (5 / 3) * r*r

	factor2 = np.exp(-np.sqrt(5) * r)

	matern_kernel = (np.log(t0) * np.log(t0))* factor1* factor2

	return matern_kernel

def exponentialKernel(xKernel, yKernel, indexI, indexJ, Dimension, k_param):
 
	#This is the Squared Exponential Kernel
	#https://nl.mathworks.com/help/stats/kernel-covariance-function-options.html
	xI = np.zeros(Dimension)
	xJ = np.zeros(Dimension)

	for i in range (0,Dimension):
		xI[i] = xKernel[indexI][i]
		xJ[i] = yKernel[indexJ][i]

	temp = subVector(xI, xJ, Dimension)

	r = 0
	for i in range(0,Dimension): 
		x1 = temp[i]
		t1 = k_param[i]
		r += (x1 * x1) /(np.log(t1) * np.log(t1))

	r = np.sqrt(r)
	t0 = k_param[Dimension]

	exponential_kernel = (np.log(t0) * np.log(t0)) * np.exp(-((r*r) / 2));

	return exponential_kernel
 
def maternKernel32(xKernel, yKernel, indexI, indexJ, Dimension, k_param):
	#This is the Matern 5/2 Kernel
	#https://nl.mathworks.com/help/stats/kernel-covariance-function-options.html
	xI = np.zeros(Dimension)
	xJ = np.zeros(Dimension)

	for i in range (0,Dimension):
		xI[i] = xKernel[indexI][i]
		xJ[i] = yKernel[indexJ][i]

	temp = subVector(xI, xJ, Dimension)

	r = 0
	for i in range(0,Dimension): 
		x1 = temp[i]
		t1 = k_param[i]
		r += (x1 * x1) /(np.log(t1) * np.log(t1))

	r = np.sqrt(r)
	t0 = k_param[Dimension]

	factor1 = 1 + np.sqrt(3) * r
	factor2 = np.exp(-np.sqrt(3) * r)

	matern_kernel = (np.log(t0) * np.log(t0)) * factor1 * factor2

	return matern_kernel

def diffMaternKernel52(xKernel, yKernel, indexI, indexJ, Dimension, k_param, index):
	xI = np.zeros(Dimension)
	xJ = np.zeros(Dimension)

	for i in range (0,Dimension):
		xI[i] = xKernel[indexI][i]
		xJ[i] = yKernel[indexJ][i]

	temp = subVector(xI, xJ, Dimension)
	
	sumatemp = 0.0
	for i in range(0, Dimension):
		sumatemp=sumatemp+temp[i]

	r = 0;
	
	for i in range(0,Dimension):
		x1 = temp[i]
		t1 = k_param[i]
		r += (x1 * x1) / (np.log(t1) * np.log(t1))
	
	r = np.sqrt(r)
	
	factor = 1
	t0 = k_param[Dimension]
	
	matern_kernel=0
	
	if index==Dimension : 
		f1 = np.sqrt(5) * r
		matern_kernel = 2 * np.log(t0)/t0 * np.exp(-f1) * ((5 / 3) * (r * r) + f1 + 1)
		
	else:
		t1 = k_param[index]
		
		x1 = temp[index]
		
		f3 = r
		
		f2 = np.exp(-np.sqrt(5) * r)
		
		f1 = t1 * np.log(t1) * np.log(t1) * np.log(t1) * r
		
		numerator = np.sqrt(5) * (np.log(t0) * np.log(t0)) * (x1 * x1) * f2 * ((5 / 3) * (r * r) + np.sqrt(5) * f3 + 1)
		
		denominator = f1
		
		matern_kernel = numerator / denominator - (np.log(t0) * np.log(t0)) * f2 * ((10 * x1 * x1) / (3 * t1 * np.log(t1) * np.log(t1) * np.log(t1)) + (np.sqrt(5) * x1 * x1) / (f1));

	if math.isnan(matern_kernel):
		matern_kernel=0

	return matern_kernel
