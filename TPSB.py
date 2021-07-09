import numpy as np

try:
    exec(open(fname).read())
except:
    exec(open('in1.py').read())

def add_Kel_to_K(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*2,inode1*2,inode2*2,inode2*2]
    lJ = [inode1*2,inode2*2,inode1*2,inode2*2]
    lN = [0,0,2,2]
    lM = [0,2,0,2]
    for I,J,N,M in zip(lI,lJ,lN,lM):
        for i in [0,1]:
            for j in [0,1]:
                K[I+i][J+j] = K[I+i][J+j] + element.Kel[N+i][M+j]

def add_Mel_to_M(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*2,inode1*2,inode2*2,inode2*2]
    lJ = [inode1*2,inode2*2,inode1*2,inode2*2]
    lN = [0,0,2,2]
    lK = [0,2,0,2]
    for I,J,N,K in zip(lI,lJ,lN,lK):
        for i in [0,1]:
            for j in [0,1]:
                M[I+i][J+j] = M[I+i][J+j] + element.Mel[N+i][K+j]

def add_Rel_to_R(element):
    inode1 = element.node1.i
    inode2 = element.node2.i
    lI = [inode1*2,inode2*2]
    lN = [0,2]
    for I,N in zip(lI,lN):
        for i in [0,1]:
            R[I+i][0] = R[I+i][0] + element.Rel[N+i][0]

nnodes = len(nodes)
dimK = nnodes*2 # caso bidimensional

elimDOFs.sort(reverse=True)
nodeMap = {}
imap = 0
for idof in range(dimK):
    if idof in elimDOFs:
        nodeMap[idof] = None
    else:
        nodeMap[idof] = imap
        imap = imap + 1
#print(nodeMap)

# Matriz de rigidez global
K = np.zeros((dimK,dimK))
for i,el in enumerate(els):
    el.update_Kloc()
    el.update_Kel()
    add_Kel_to_K(el)
for DOF in elimDOFs:
    K = np.delete(K, DOF, 0)
    K = np.delete(K, DOF, 1)

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

# calcula R inicial
R = np.zeros((dimK,1))
for i,el in enumerate(els):
    el.update_Rel()
    add_Rel_to_R(el)
for DOF in elimDOFs:
    R = np.delete(R, DOF, 0)

# Forca periodica
w = 2.0*np.pi/T
t = np.linspace(0.0, tot_t, nsteps+1)
FF = np.vstack([Fi*np.sin(w*t) for Fi in F])

##LL = np.array([0])
##ee = np.array([0])
##SS = np.array([0])
LL = np.zeros((nsteps+1,1))
ee = np.zeros((nsteps+1,1))
SS = np.zeros((nsteps+1,1))

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

def dprint(*args):
    print(*args)

DEBUG = False
if DEBUG:
    nsteps = 251
    progress_bar = lambda *args: None
else:
    dprint = lambda *args: None

step = 0
iters = []
STBK = False
while step < nsteps:

    step = step + 1
    
    ia = step - 1
    U[step] = U[ia]
    
    niter = 0
    Rnorm = np.inf
    
    for el in els:
        E = el.material.E
        if isinstance(E, Plasticity1D):
            E.backup()
    
    while Rnorm > 0.001 and niter < 500:
    
        niter = niter + 1
        
        auxM1 = ic[0] * (U[ia] - U[step]) + ic[2] * V[ia] + ic[3] * A[ia]
        ##auxM1 = ic[0] * (U[step] - U[ia]) + ic[2] * V[ia] + ic[3] * A[ia]
        auxM2 = np.matmul(M,auxM1)
        
        ##auxC1 = ic[1] * U[ia] + ic[4] * V[ia] + ic[5] * A[ia]
        ##auxC2 = np.matmul(C,auxC1)
        
        FI = auxM2 ##+ auxC2
        
        Fext = np.transpose([FF[:,step]])
        FE = Fext + FI
        
        dR = FE - R
        if niter == 1:
            dprint()
            dprint("### STEP {} ###".format(step))
        dprint()
        dprint("ITER {}".format(niter))
        dprint()
        dprint("- Forca efetiva")
        dprint(dR,FE,R)
        dprint()
        
        for i,el in enumerate(els):
            E = el.material.E
            if isinstance(E, Plasticity1D):
                E.restore()
            el.update_Kloc()
            if isinstance(E, Plasticity1D):
                
                if STBK:
            
                    dprint("----- STEPBACK")
                    
                    if niter == 1:
                        el.e_i = E.e
                        el.S_i = E.S
                        dprint("DEF. INICIAL = {}".format(E.e))
                        dprint("TENSAO INICIAL = {}".format(E.S))
                        dprint("FORCA INICIAL = {}".format(E.S*el.material.A))
                        
                    inode1x = 2*el.node1.i
                    inode1y = 2*el.node1.i+1
                    Floc1x = FE[nodeMap[inode1x]] if isinstance(nodeMap[inode1x],int) else 0
                    Floc1y = FE[nodeMap[inode1y]] if isinstance(nodeMap[inode1y],int) else 0
                    Floc1 = Floc1x*el.cos + Floc1y*el.sen
                    inode2x = 2*el.node2.i
                    inode2y = 2*el.node2.i+1
                    Floc2x = FE[nodeMap[inode2x]] if isinstance(nodeMap[inode2x],int) else 0
                    Floc2y = FE[nodeMap[inode2y]] if isinstance(nodeMap[inode2y],int) else 0
                    Floc2 = Floc2x*el.cos + Floc2y*el.sen
                    Floc = Floc2 - Floc1
                    
                    dprint("Forcas locais no 1")
                    dprint(Floc1x,Floc1y,Floc1)
                    dprint("Forcas locais no 2")
                    dprint(Floc2x,Floc2y,Floc2)
                    dprint("Forca local total")
                    dprint(Floc)
                    
                    S_f = Floc/el.material.A
                    e_i = el.e_i
                    S_i = el.S_i
                    
                    ##if DEBUG:
                    ##    e_f = S_f/E1
                    ##    if S_f > Sy:
                    ##        e_f = Sy/E1 + (S_f-Sy)/E2
                    ##else:
                    ##    E.backup()
                    ##    E.update(S = S_f)
                    ##    e_f = E.e
                    ##    E.restore()
                    E.backup()
                    E.update(S = S_f)
                    e_f = E.e
                    E.restore()
                    
                    de = e_f - e_i
                    dS = S_f - S_i
                    Esec = abs(dS/de)
                    Esec = np.where(np.isnan(Esec),0,Esec)
                    
                    dprint("NOVA TENSAO = {}".format(S_f))
                    dprint("NOVA DEF. = {}".format(e_f))
                    dprint(dS,de,Esec)
                    
                    Ksec = Esec*el.material.A/el.l0
                    el.Ksec = Ksec[0]
                    dprint(el.Ksec,el.Kloc)
                    el.Kloc = np.where(el.Ksec == 0,el.Kloc,el.Ksec)
                    dprint(el.Ksec,el.Kloc)
                
                ##if niter == 1:
                ##    ##LL = np.append(LL,(el2.length-el2.l0)/el2.l0)
                ##    ee = np.append(ee,E.e)
                ##    SS = np.append(SS,E.S)
                ee[step] = E.e
                SS[step] = E.S
                
            dprint(el.Kloc)
        
        # atualiza K
        K = np.zeros((dimK,dimK))
        for i,el in enumerate(els):
            el.update_Kel()
            add_Kel_to_K(el)
        for DOF in elimDOFs:
            K = np.delete(K, DOF, 0)
            K = np.delete(K, DOF, 1)
        
        KE = K + ic[0] * M ##+ ic[1] * C
        
        # incrementa deslocamentos
        KEinv = np.linalg.inv(KE)
        dU = np.matmul(KEinv,dR)
        U[step] = U[step] + dU
        
        #dprint()
        #dprint("----- NEWTON RAPHSON")
        #dprint("Matriz de rigidez")
        #dprint(KE)
        #dprint("Inversa")
        #dprint(KEinv)
        #dprint("Deslocamentos")
        #dprint(dU)
        
        # atualiza nos a partir de U
        for node in nodes:
            DOFx = nodeMap[2*node.i]
            DOFy = nodeMap[2*node.i+1]
            if DOFx != None:
                node.x = node.x0 + U[step][DOFx][0]
            if DOFy != None:
                node.y = node.y0 + U[step][DOFy][0]
        
        ##dprint("- Coordenadas dos nos:")
        ##for inode, node in enumerate(nodes):
        ##    dprint("No {0}: ({1:g}; {2:g})".format(inode,node.x,node.y))
        ##dprint()
        
        # atualiza R
        R = np.zeros((dimK,1))
        for i,el in enumerate(els):
            el.update_Rel()
            add_Rel_to_R(el)
        for DOF in elimDOFs:
            R = np.delete(R, DOF, 0)
        
        # reavalia residuo
        dR = FE - R
        Rnorm = np.linalg.norm(dR)
        #dprint("- Norma do residuo: {:g}".format(Rnorm))
        #dprint()
    
    # Atualiza aceleracoes e velocidades
    A[step] = ic[0] * (U[step] - U[ia]) - ic[2] * V[ia] - ic[3] * A[ia]
    V[step] = V[ia] + ic[6] * A[ia] + ic[7] * A[step]
    
    ##LL = np.append(LL,(el2.length-el2.l0)/el2.l0)
    ##ee = np.append(ee,el2.material.E.e)
    ##SS = np.append(SS,el2.material.E.S)
    
    ##LL = np.append(LL,(el1.length-el1.l0)/el1.l0)
    ##ee = np.append(ee,el1.material.E.e)
    ##SS = np.append(SS,el1.material.E.S)
    
    iters.append(niter)
    progress_bar(step,nsteps)

E = lambda i: (SS[i+1]-SS[i])/(ee[i+1]-ee[i])

import matplotlib.pyplot as plt

#if STBK or DEBUG:
if True:

    fig = plt.figure(figsize=(10,5))
    plt.plot(ee,SS,'-o')
    
    Eaux = Plasticity1D(E1,E2,Sy)
    e_ana = np.hstack([
        np.linspace(0,ee.max(),10000),
        np.linspace(ee.max(),ee.min(),10000),
        np.linspace(ee.min(),0,10000),
    ])
    #e_ana = np.insert(e_ana,len(e_ana)-sum(e_ana>ecrit),ecrit)
    #S_ana = np.array([S_x_e(e) for e in e_ana])
    S_ana = np.array([])
    for e in e_ana:
        Eaux.update(e = e)
        S_ana = np.append(S_ana, Eaux.S)
    plt.plot(e_ana,S_ana,'-')
    
    plt.xlabel("e")
    plt.ylabel("S")
    plt.show()

if DEBUG == False:

    data = np.transpose(U)
    t = np.linspace(0, nsteps*Dt, nsteps+1)
    
    N = int(1.0*nsteps)
    ###N = int(1.0*nsteps-510)
    
    fig = plt.figure(figsize=(10,5))
    plt.subplot(311)
    plt.plot(t[:N],FF[1][:N],'-')
    plt.xlabel("t")
    plt.ylabel("F")
    plt.subplot(312)
    plt.plot(t[:N],data[0][0][:N],'-')
    plt.xlabel("t")
    plt.ylabel("x(1)")
    plt.subplot(313)
    plt.plot(t[:N],data[0][1][:N],'-')
    plt.xlabel("t")
    plt.ylabel("x(2)")
    ##plt.subplot(514)
    ##plt.plot(t[:N],ee[:N],'-')
    ##plt.xlabel("t")
    ##plt.ylabel("e")
    ##plt.subplot(515)
    ##plt.plot(ee[:N],LL[:N],'-')
    ##plt.xlabel("e")
    ##plt.ylabel("S")
    plt.show()
    print()
