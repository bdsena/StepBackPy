import numpy as np
import matplotlib.pyplot as plt

data1 = np.load('out.npz')
data2 = np.load('outSB.npz')
SS1 = data1['SS'][1]#[-4:]
ee1 = data1['ee'][1]#[-4:]
SS2 = data2['SS'][1]#[-4:]
ee2 = data2['ee'][1]#[-4:]

fig = plt.figure()
plt.plot(ee1,SS1,'-o')
plt.plot(ee2,SS2,'-o')
plt.xlabel("e")
plt.ylabel("S")
plt.show()
