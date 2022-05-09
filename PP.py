import sys
import numpy as np
from copy import deepcopy

try:
    exec(open(sys.argv[1]).read())
except:
    exec(open('in2.py').read())

def add_Kel_to_K(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*3,inode1*3,inode2*3,inode2*3]
    lJ = [inode1*3,inode2*3,inode1*3,inode2*3]
    lN = [0,0,3,3]
    lM = [0,3,0,3]
    for I,J,N,M in zip(lI,lJ,lN,lM):
        for i in [0,1,2]:
            for j in [0,1,2]:
                K[I+i][J+j] = K[I+i][J+j] + element.Kel[N+i][M+j]

def add_Mel_to_M(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*3,inode1*3,inode2*3,inode2*3]
    lJ = [inode1*3,inode2*3,inode1*3,inode2*3]
    lN = [0,0,3,3]
    lK = [0,3,0,3]
    for I,J,N,K in zip(lI,lJ,lN,lK):
        for i in [0,1,2]:
            for j in [0,1,2]:
                M[I+i][J+j] = M[I+i][J+j] + element.Mel[N+i][K+j]

def add_Rel_to_R(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*3,inode2*3]
    lN = [0,3]
    for I,N in zip(lI,lN):
        for i in [0,1,2]:
            R[I+i][0] = R[I+i][0] + element.Rel[N+i][0]

nnodes = len(nodes)
dimK = nnodes*3 # caso bidimensional

elimDOFs.sort(reverse=True)
nodeMap = {}
imap = 0
for idof in range(dimK):
    if idof in elimDOFs:
        nodeMap[idof] = None
    else:
        nodeMap[idof] = imap
        imap = imap + 1

# Matriz de massa global
M = np.zeros((dimK,dimK))
for i,el in enumerate(els):
    add_Mel_to_M(el)
for DOF in elimDOFs:
    M = np.delete(M, DOF, 0)
    M = np.delete(M, DOF, 1)

# Constantes de integracao
alpha = 0.25
beta = 0.5
ic = np.zeros(8)
ic[0] = 1.0 / (alpha * Dt * Dt)
ic[1] = beta / (alpha * Dt)
ic[2] = 1.0 / (alpha * Dt)
ic[3] = 1.0 / (2.0 * alpha) - 1.0
ic[4] = beta / alpha - 1.0
ic[5] = Dt * (beta / alpha - 2.0) / 2.0
ic[6] = Dt * (1.0 - beta)
ic[7] = beta * Dt

nsteps = int(tot_t / Dt)
U = np.zeros((nsteps+1,len(F),1))
V = np.zeros((nsteps+1,len(F),1))
A = np.zeros((nsteps+1,len(F),1))

# Forca periodica
w = 2.0*np.pi/T
t = np.linspace(0.0, tot_t, nsteps+1)
#FF = np.vstack([Fi*np.sin(w*t) for Fi in F])
#ramp_func = np.linspace(0.0,1.0,nsteps+1)
ramp_func = np.hstack([
    np.linspace( 0.00, 1.00,int(nsteps/4)+1),
    np.linspace( 1.00,-1.05,int(nsteps/2)+1)[1:],
    np.linspace(-1.05, 0.05,int(nsteps/4)+1)[1:],
])
FF = np.vstack([Fi*ramp_func for Fi in F])
print(FF[-1])

ee = np.zeros((len(els),1))
SS = np.zeros((len(els),1))

import time
import sys

t0 = time.time()
def progress_bar(i,N):
    bar_len = 60
    filled_len = int(round(bar_len*i/N))
    percents = round(100*i/N,1)
    bar = '='*filled_len + '-'*(bar_len-filled_len)
    t = round(time.time()-t0,2)
    sys.stdout.write('[%s] %s%s (%s/%s) Elapsed time: %ss\r' % (bar,percents,'%',i,N,t))
    sys.stdout.flush()

class Iters:
    def __init__(self):
        self.values = []
    def __repr__(self):
        return repr(self.values)
    def add(self, new):
        self.values.append(new)
    def sum(self):
        return sum([sum(i) for i in self.values])
    def step(self):
        return [sum(i) for i in self.values]
    def stbk(self):
        return [len(i)-1 for i in self.values]
    
R = np.zeros((dimK,1))
for el in els:
    el.set_STBK(STBK)

step = 0
maxiter = 200
maxloop = 5
iters = Iters()
while step < nsteps:

    step = step + 1
    ia = step - 1
    
    for el in els:
        el.backup()
    
    step_iters = []
    
    iloop = 0
    NRP_Loop = True
    while NRP_Loop and iloop < maxloop:
    
        U[step] = U[ia]
        
        niter = 0
        Rnorm = np.inf

        iloop = iloop + 1
        Rbak = np.copy(R)
        
        while Rnorm > 1e-4 and niter < maxiter:
        
            niter = niter + 1
            #print('============================= iter ',niter)

            # atualiza R
            R = np.zeros((dimK,1))
            for i,el in enumerate(els):
                el.update_length()
                el.update_Rloc()
                el.update_Rel()
                add_Rel_to_R(el)
            for DOF in elimDOFs:
                R = np.delete(R, DOF, 0)
            
            # atualiza K
            K = np.zeros((dimK,dimK))
            for i,el in enumerate(els):
                el.update_Kloc(step)
                el.update_Kel()
                add_Kel_to_K(el)
            for DOF in elimDOFs:
                K = np.delete(K, DOF, 0)
                K = np.delete(K, DOF, 1)
            
            auxM1 = ic[0] * (U[ia] - U[step]) + ic[2] * V[ia] + ic[3] * A[ia]
            auxM2 = np.matmul(M,auxM1)
            FI = auxM2
            
            Fext = np.transpose([FF[:,step]])
            FE = Fext + FI
            
            dR = FE - R
            ##print('# pre')
            ##print('FE')
            ##print(FE)
            ##print('R')
            ##print(R)
            ##print('dR')
            ##print(dR)
            
            KE = K + ic[0] * M
            #print('KE')
            #print(KE)
            
            # incrementa deslocamentos
            KEinv = np.linalg.inv(KE)
            #print('KEinv')
            #print(KEinv)
            dU = np.matmul(KEinv,dR)
            #print('dR')
            #print(dR)
            #print('dU')
            #print(dU)
            U[step] = U[step] + dU
            
            # atualiza nos a partir de U
            for node in nodes:
                DOFx = nodeMap[3*node.i]
                DOFy = nodeMap[3*node.i+1]
                DOFa = nodeMap[3*node.i+2]
                if DOFx != None:
                    node.x = node.x0 + U[step][DOFx][0]
                if DOFy != None:
                    node.y = node.y0 + U[step][DOFy][0]
                if DOFa != None:
                    node.a = node.a0 + U[step][DOFa][0]
                    node.ainc = dU[DOFa][0]
                #print('i = ',node.i)
                #print('x = ',node.x)
                #print('y = ',node.y)
                #print('a = ',node.a)
                #print('ainc = ',node.ainc)
            
            # reavalia residuo
            Rnorm = np.linalg.norm(dU)
        
        step_iters.append(niter)
        
        NRP_Loop = False
        for el in els:
            NRP_Loop_el = el.NRP_check()
            NRP_Loop = NRP_Loop or NRP_Loop_el
        
        if NRP_Loop:
            
            R = np.copy(Rbak)
            
            # restaura nos a partir de U
            for node in nodes:
                DOFx = nodeMap[3*node.i]
                DOFy = nodeMap[3*node.i+1]
                DOFa = nodeMap[3*node.i+2]
                if DOFx != None:
                    node.x = node.x0 + U[ia][DOFx][0]
                if DOFy != None:
                    node.y = node.y0 + U[ia][DOFy][0]
                if DOFa != None:
                    node.a = node.a0 + U[ia][DOFa][0]
            
            for el in els:
                el.restore()
        
    if iloop == maxloop:
        print()
        print("CONV ERROR STEP {}".format(step))
        step = nsteps
    
    # Atualiza aceleracoes e velocidades
    A[step] = ic[0] * (U[step] - U[ia]) - ic[2] * V[ia] - ic[3] * A[ia]
    V[step] = V[ia] + ic[6] * A[ia] + ic[7] * A[step]
    
    new_ee = []
    new_SS = []
    for i,el in enumerate(els):
        new_e = el.k1
        new_S = -el.Rloc[2][0]
        new_ee.append(new_e)
        new_SS.append(new_S)
    ee = np.hstack([ee, np.transpose([new_ee])])
    SS = np.hstack([SS, np.transpose([new_SS])])
    
    iters.add(step_iters)
    progress_bar(step,nsteps)

print()
print(iters.values)
print(iters.sum())

desl = np.transpose(U)[0]
t = np.linspace(0, nsteps*Dt, nsteps+1)

xnodes = np.ones((nsteps+1, nnodes))
ynodes = np.ones((nsteps+1, nnodes))
anodes = np.ones((nsteps+1, nnodes))
for i,node in enumerate(nodes):
    idofx = nodeMap[3*i]
    if idofx != None:
        xnodes[:,i] = xnodes[:,i]*node.x0 + U[:,idofx,0]
    else:
        xnodes[:,i] = xnodes[:,i]*node.x0
    idofy = nodeMap[3*i+1]
    if idofy != None:
        ynodes[:,i] = ynodes[:,i]*node.y0 + U[:,idofy,0]
    else:
        ynodes[:,i] = ynodes[:,i]*node.y0
    idofa = nodeMap[3*i+2]
    if idofa != None:
        anodes[:,i] = anodes[:,i]*node.a0 + U[:,idofa,0]
    else:
        anodes[:,i] = anodes[:,i]*node.a0

dofs = np.array([nodeMap[k] for k in nodeMap], dtype=float)
np.savez('outPP.npz',t=t,desl=desl,FF=FF,ee=ee,SS=SS,xnodes=xnodes,ynodes=ynodes,anodes=anodes,dofs=dofs,ndof=3)
