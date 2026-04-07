import numpy as np
import cma
from testFunctions import * 
from mathFunctions import * 
import sys
from ctypes import *
from scipy.stats import qmc
import time
from pandas import read_csv
import pipeline0

def cmaesCode():
	Dimension=6
	x = 6 * [1]  # initial solution
	print(x)
	sigma0 = 0.1	# initial standard deviation to sample new solutions
	xopt, es = cma.fmin2(hartmann6D, x, sigma0, options=None, args=[Dimension])
	
	return


def mainFunction(x, Dimension, run_id=0,
                 base_dir="base", run_prefix="run",
                 size_x=20, size_y=20, size_z=20,
                 num_modes=10, energy_range=4, exhaustiveness=16):

    # x is assumed in [0,1]^Dimension for the first 3 dims
    print("x values", x)

    xMod = np.zeros(Dimension)

    # Correct receptor bounds (span and min) from your structure
    # span_x=61.072, span_y=56.715, span_z=103.493
    # min_x=-5.756, min_y=14.686, min_z=24.578
    xMod[0] = (61.072) * x[0] - 5.756
    xMod[1] = (56.715) * x[1] + 14.686
    xMod[2] = (103.493) * x[2] + 24.578

    print("xMod values", xMod)

    # Treat xMod[0:3] as the docking box center
    center_x = float(xMod[0])
    center_y = float(xMod[1])
    center_z = float(xMod[2])

    # Run docking pipeline and return the mean affinity
    mean_val = pipeline0.mainPipeline(
        base_dir=base_dir,
        run_id=run_id,
        run_prefix=run_prefix,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        num_modes=num_modes,
        energy_range=energy_range,
        exhaustiveness=exhaustiveness
    )

    print("Returned mean affinity", mean_val)
    return mean_val


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
		ytrain[i]=mainFunction(xtrain[i], Dimension,i)

	return xtrain, ytrain

#UCBParallel1KernelScalable
def mainTest():
	# Start timer
	start_time = time.perf_counter()

	Dimension = 3
	n_points_Start = 5 #inital points
	evalsMethod = 300000 #30,000 75,000 150,000 300,000
	method = 1 #DIRECT,CMAES in the selection of the point
	iterationsMax = 100 #100,150,300,500
	subset = 50 #subset of points 
	threadsParallel = 5
	typeKernel=0#Matern 52, Exponential, Matern 32
	runName=1#variable used to write the results
	noiseLevel = 1e-8

	tries=threadsParallel
	
	n_points = n_points_Start
	
	xtrain, ytrain = generateInitialPoints(n_points, Dimension, n_points_Start)
	#xtrain, ytrain, n_points=loadDataset()
	
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
		nextPoint=mainFunction(xopt, Dimension, n_points)
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
	
	
	
			