import numpy as np
from testFunctions import * 
from mathFunctions import *
from pandas import read_csv 
import pandas as pd
from scipy.stats import qmc

folderName="data/"

def mainFunction(x, Dimension):
	resultTemp=sphereFunctionProblemEvaluation(x, Dimension)
	print("x values", x)
	xMod=np.zeros(Dimension)
	xMod[0]=(60.903)*x[0]-5.756
	xMod[1]=(54.731)*x[1]+16.108
	xMod[2]=(106.062)*x[2]+21.029
	print("xMod values", xMod)
	#xMod[0]=round(x[0] * ((0.9-0.1)/0.05)) *0.05+0.1
	#xMod[1]=round(x[1] * ((120-5)/5)) *5+5
	#xMod[2]=round(x[2] * ((25-0)/0.5)) *0.5+0
	#print("xMod values", xMod)
	print("possible Result %.4f"%(resultTemp))
	result=float(input("Press Value + Enter to continue... "))
	return result

def generateInitialPoints(n_points, Dimension):

	sampler = qmc.LatinHypercube(d=Dimension)
	sample = sampler.random(n=n_points)
	print(sample)

	#xtrain = np.zeros((n_points, Dimension))
	xtrain = np.copy(sample)
	pd.DataFrame(xtrain).to_csv(folderName+"data.csv", header=None, index =None)
	ytrain = np.zeros(n_points)
	for i in range(0,n_points):
		#for j in range (0, Dimension):
		#	xtrain[i][j]=np.random.rand()
		ytrain[i]=mainFunction(xtrain[i], Dimension)
	
	pd.DataFrame(ytrain).to_csv(folderName+"labels.csv", header=None, index =None)
	
	
	return
	
if __name__ == "__main__":
	generateInitialPoints(5, 3)