import numpy as np
x=np.random.rand(3)
print("x values", x)
xMod=np.zeros(3)
xMod[0]=(0.9-0.1)*x[0]+0.1
xMod[1]=(120-5)*x[1]+5
xMod[2]=(25-0)*x[2]+0
print("xMod values", xMod)
xMod[0]=round(x[0] * ((0.9-0.1)/0.05)) *0.05+0.1
xMod[1]=round(x[1] * ((120-5)/5)) *5+5
xMod[2]=round(x[2] * ((25-0)/0.5)) *0.5+0
print("xMod values", xMod)


	