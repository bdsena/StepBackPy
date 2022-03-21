import sys
import numpy as np
import matplotlib.pyplot as plt

data1 = np.load('out.npz')
data2 = np.load('outSB.npz')

if 'desl' in sys.argv:
    inode1 = 2
    inode2 = 2
    idof1 = 0
    idof2 = 0
    dofs1 = data1['dofs']
    ndof1 = data1['ndof']
    dofs2 = data2['dofs']
    ndof2 = data2['ndof']
    idesl1 = int(dofs1[inode1*ndof1+idof1])
    idesl2 = int(dofs2[inode2*ndof2+idof2])
    t1 = data1['t']
    desl1 = data1['desl'][idesl1]
    t2 = data2['t']
    desl2 = data2['desl'][idesl2]
    
    fig = plt.figure()
    plt.plot(t1,desl1)
    plt.plot(t2,desl2)
    plt.xlabel("t")
    plt.ylabel("x({0})".format(idof1))
    plt.show()

else:
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
