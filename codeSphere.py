import numpy as np
import cma
from scipy.optimize import minimize
from scipy.optimize import Bounds
from testFunctions import * 
from mathFunctions import * 
from gaussian import * 
import sys
from ctypes import *
# Load DLL into memory.
from scipy.stats import qmc
import time
from pandas import read_csv

def cmaesCode():
	Dimension=6
	x = 6 * [1]  # initial solution
	print(x)
	sigma0 = 0.1	# initial standard deviation to sample new solutions
	xopt, es = cma.fmin2(hartmann6D, x, sigma0, options=None, args=[Dimension])
	
	return

def mainFunction(x, Dimension):
	resultTemp=sphereFunctionProblemEvaluation(x, Dimension)
	print("x values", x)
	xMod=np.zeros(Dimension)
	xMod[0]=(0.9-0.1)*x[0]+0.1
	xMod[1]=(120-5)*x[1]+5
	xMod[2]=(25-0)*x[2]+0
	print("xMod values", xMod)
	xMod[0]=round(x[0] * ((0.9-0.1)/0.05)) *0.05+0.1
	xMod[1]=round(x[1] * ((120-5)/5)) *5+5
	xMod[2]=round(x[2] * ((25-0)/0.5)) *0.5+0
	print("xMod values", xMod)
	print("possible Result %.4f"%(resultTemp))
	result=float(input("Press Enter to continue... "))
	return result

def loadDataset() :
	
	# data used for the predictions
	dfData = read_csv("./data/data.csv", header=None, sep=',', dtype=float)
	dfLabels = read_csv("./data/labels.csv", header=None, dtype=float)
	n_points=len(dfLabels.values.ravel())
	
	return dfData.values, dfLabels.values.ravel(), n_points
	
def generateInitialPoints(n_points, Dimension, n_points_Start):

	sampler = qmc.LatinHypercube(d=Dimension)
	sample = sampler.random(n=n_points_Start)
	print(sample)

	#xtrain = np.zeros((n_points, Dimension))
	xtrain = np.copy(sample)
	ytrain = np.zeros(n_points)
	for i in range(0,n_points):
		#for j in range (0, Dimension):
		#	xtrain[i][j]=np.random.rand()
		ytrain[i]=mainFunction(xtrain[i], Dimension)

	return xtrain, ytrain

#UCBParallel1KernelScalable
def mainTest():
	# Start timer
	start_time = time.perf_counter()

	Dimension = 3
	n_points_Start = 3 #inital points
	evalsMethod = 50000 #30,000 75,000 150,000 300,000
	method = 1 #DIRECT,CMAES in the selection of the point
	iterationsMax = 500 #100,150,300,500
	subset = 50 #subset of points 
	threadsParallel = 5
	typeKernel=0#Matern 52, Exponential, Matern 32
	runName=1#variable used to write the results
	noiseLevel = 1e-8

	tries=threadsParallel
	
	n_points = n_points_Start
	
	#xtrain, ytrain = generateInitialPoints(n_points, Dimension, n_points_Start)
	xtrain, ytrain, n_points=loadDataset()
	
	bestValue = 1e10

	np.set_printoptions(precision=6)
	
	nameFile=str(runName)+'_output.txt'
	
	with open(nameFile, 'w') as f:
		f.close()
	
	for n in range(0,n_points):
			if ytrain[n] <bestValue:
				bestValue = ytrain[n]
	
	with open(nameFile, 'a') as f:
		for i in range(0, n_points):
			f.write("%0.6f\t%0.6f\t"%(bestValue,ytrain[i]))
			for j in range(0,Dimension):
				if j == Dimension-1:
					f.write("%0.6f\n"%(xtrain[i][j]))
				else:
					f.write("%0.6f\t"%(xtrain[i][j]))
		f.close()
	hllDll = CDLL (".\Solver.dll")
	hllDll.getNextPoint.restype=POINTER(c_double* Dimension)
	
	for iterations in range(0,iterationsMax):
		
		for n in range(0,n_points):
			if ytrain[n] <bestValue:
				bestValue = ytrain[n]
				bestPoint = xtrain[n]
		
		print("bestValue %.6f"%bestValue)
		
		#hllDll.showDouble(c_double(10.0))
		
		yPoints=ytrain.ctypes.data_as(POINTER(c_double))
		#hllDll.showVector(yPoints, c_int(n_points))
		
		DOUBLE = c_double
		PDOUBLE = POINTER(DOUBLE)
		PPDOUBLE = POINTER(PDOUBLE)
		INT = c_int
		
		# An array of doubles can be passed to a function that takes double*.
		DBL5ARR = DOUBLE * Dimension
		# An array of double* can be passed to your function as double**.
		PDBL4ARR = PDOUBLE * n_points

		# Declare double* array.
		xPoints = PDBL4ARR()
		for i in range(n_points):
		# fill out each pointer with an array of doubles.
			xPoints[i] = DBL5ARR()
			for j in range(Dimension):
				xPoints[i][j] = c_double(xtrain[i][j])
		
		#hllDll.showMatrix(xPoints, c_int(n_points), c_int(Dimension))
		
# (int n_points, int Dimension, int threadsParallel, double** xtrain,
	# double* ytrain, int subset, int tries, int method, int typeKernel, double bestValue,
	# int evalsMethod, int runName, double noiseLevel);
		
		
		result=hllDll.getNextPoint(c_int(n_points), c_int(Dimension), c_int(threadsParallel), xPoints,
			yPoints, c_int(subset), c_int(tries), c_int(method), c_int(typeKernel), c_double(bestValue), 
			c_int(evalsMethod), c_int(runName), c_double(noiseLevel), c_int(iterations)).contents
		
		xopt=np.zeros(Dimension)
		for i in range(Dimension):
			xopt[i]=result[i]
					
		print("xopt", xopt)
		nextPoint=mainFunction(xopt, Dimension)
		#print(xtrain)
		#print(ytrain)
		
		xtrain=np.append(xtrain, [xopt], axis=0)
		ytrain=np.append(ytrain, nextPoint)

		#print(xtrain)
		#print(ytrain)

		
		
		with open(nameFile, 'a') as f:
			f.write("%0.6f\t%0.6f\t"%(bestValue, nextPoint))
			for i in range(0,Dimension):
				if i == Dimension-1:
					f.write("%0.6f\n"%(xtrain[n_points][i]))
				else:
					f.write("%0.6f\t"%(xtrain[n_points][i]))
			f.close()
		
		n_points+=1
		print(n_points)
	
	# End timer
	end_time = time.perf_counter()
	
	# Calculate elapsed time
	elapsed_time = end_time - start_time
	
	print("Elapsed time: ", elapsed_time)
	
	timeFile=str(runName)+"_outputTime.txt" 
	with open(timeFile, 'w') as f:
		f.write("%0.6f\n"%(elapsed_time))
		f.close()
	return

def severalRuns():
	
	Dimension = 20
	n_points_Start = 10 #inital points
	evalsMethod = 50000 #30,000 75,000 150,000 300,000
	method = 1 #DIRECT,CMAES in the selection of the point
	iterationsMax = 500 #100,150,300,500
	subset = 50 #subset of points 
	threadsParallel = 10
	typeKernel=0#Matern 52, Exponential, Matern 32
	runName=1#variable used to write the results
	noiseLevel = 1e-8
	
	for run in range(5,20):
			# Start timer
		start_time = time.perf_counter()

		runName=run

		tries=threadsParallel
		
		n_points = n_points_Start
		
		xtrain, ytrain = generateInitialPoints(n_points, Dimension, n_points_Start)
		
		bestValue = 1e10

		np.set_printoptions(precision=6)
		
		nameFile=str(runName)+'_output.txt'
		
		with open(nameFile, 'w') as f:
			f.close()
		
		for n in range(0,n_points):
				if ytrain[n] <bestValue:
					bestValue = ytrain[n]
		
		with open(nameFile, 'a') as f:
			for i in range(0, n_points):
				f.write("%0.6f\t%0.6f\t"%(bestValue,ytrain[i]))
				for j in range(0,Dimension):
					if j == Dimension-1:
						f.write("%0.6f\n"%(xtrain[i][j]))
					else:
						f.write("%0.6f\t"%(xtrain[i][j]))
			f.close()
		hllDll = CDLL (".\Solver.dll")
		hllDll.getNextPointEI.restype=POINTER(c_double* Dimension)
		
		for iterations in range(0,iterationsMax):
			
			for n in range(0,n_points):
				if ytrain[n] <bestValue:
					bestValue = ytrain[n]
					bestPoint = xtrain[n]
			
			print("bestValue %.6f"%bestValue)
			
			#hllDll.showDouble(c_double(10.0))
			
			yPoints=ytrain.ctypes.data_as(POINTER(c_double))
			#hllDll.showVector(yPoints, c_int(n_points))
			
			DOUBLE = c_double
			PDOUBLE = POINTER(DOUBLE)
			PPDOUBLE = POINTER(PDOUBLE)
			INT = c_int
			
			# An array of doubles can be passed to a function that takes double*.
			DBL5ARR = DOUBLE * Dimension
			# An array of double* can be passed to your function as double**.
			PDBL4ARR = PDOUBLE * n_points

			# Declare double* array.
			xPoints = PDBL4ARR()
			for i in range(n_points):
			# fill out each pointer with an array of doubles.
				xPoints[i] = DBL5ARR()
				for j in range(Dimension):
					xPoints[i][j] = c_double(xtrain[i][j])
			
			#hllDll.showMatrix(xPoints, c_int(n_points), c_int(Dimension))
			
	# (int n_points, int Dimension, int threadsParallel, double** xtrain,
		# double* ytrain, int subset, int tries, int method, int typeKernel, double bestValue,
		# int evalsMethod, int runName, double noiseLevel);
			
			
			result=hllDll.getNextPointEI(c_int(n_points), c_int(Dimension), c_int(threadsParallel), xPoints,
				yPoints, c_int(subset), c_int(tries), c_int(method), c_int(typeKernel), c_double(bestValue), 
				c_int(evalsMethod), c_int(runName), c_double(noiseLevel), c_int(iterations)).contents
			
			xopt=np.zeros(Dimension)
			for i in range(Dimension):
				xopt[i]=result[i]
						
			print("xopt", xopt)
			nextPoint=mainFunction(xopt, Dimension)
			#print(xtrain)
			#print(ytrain)
			
			xtrain=np.append(xtrain, [xopt], axis=0)
			ytrain=np.append(ytrain, nextPoint)

			#print(xtrain)
			#print(ytrain)

			
			
			with open(nameFile, 'a') as f:
				f.write("%0.6f\t%0.6f\t"%(bestValue, nextPoint))
				for i in range(0,Dimension):
					if i == Dimension-1:
						f.write("%0.6f\n"%(xtrain[n_points][i]))
					else:
						f.write("%0.6f\t"%(xtrain[n_points][i]))
				f.close()
			
			n_points+=1
			print(n_points)
		
		# End timer
		end_time = time.perf_counter()
		
		# Calculate elapsed time
		elapsed_time = end_time - start_time
		
		print("Elapsed time: ", elapsed_time)
		
		timeFile=str(runName)+"_outputTime.txt" 
		with open(timeFile, 'w') as f:
			f.write("%0.6f\n"%(elapsed_time))
			f.close()
	return
if __name__ == "__main__":
	mainTest()  
	
	
	
			