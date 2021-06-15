import numpy as np

class Plasticity1D:
  
    def __init__(self, E, Et, Sy):
    
        #Dados do material
        self.E = E    # Modulo de elasticidade
        self.Et = Et  # Modulo de elasticidade tangente
        self.Sy = Sy  # Tensao de escoamento
        self.H = E * Et / (E - Et) # Modulo de encruamento
    
        #Variaveis de estado
        self.e_in = 0  # Parcela inelastica das deformacoes
        self.beta = 0  # Variavel de controle da plasticidade
        self.Ec = E    # Modulo de elasticidade atual
        self.S = 0     # Tensao atual
        self.e = 0     # Deformacao atual
        
        self.e_in_bak = 0
        self.beta_bak = 0
        self.Ec_bak = E
        self.S_bak = 0
        self.e_bak = 0

    def backup(self):
        self.e_in_bak = self.e_in
        self.beta_bak = self.beta
        self.Ec_bak = self.Ec
        self.S_bak = self.S
        self.e_bak = self.e
    def restore(self):
        self.e_in = self.e_in_bak
        self.beta = self.beta_bak
        self.Ec = self.Ec_bak
        self.S = self.S_bak
        self.e = self.e_bak

    # Atualiza o estado de tensao
    def update(self, **kwargs):
        if 'e' in kwargs:
            # Calcula como se fosse regime elastico
            self.e = kwargs.get('e')
            e_e = self.e - self.e_in
            self.S = self.E * e_e
            self.Ec = self.E
            # Verifica e faz a correcao plastica
            f = abs(self.S - self.beta) - self.Sy
            if f > 0:
                # Corrige a tensao e incrementa a parcela inelastica da deformacao
                self.Ec = self.Et
                r = self.S / abs(self.S)
                l = f / (self.E + self.H)
                self.S = self.S - self.E * l * r
                self.e_in = self.e_in + l * r
                self.beta = self.beta + l * self.H * r
        elif 'S' in kwargs:
            # Calcula como se fosse regime elastico
            self.S = kwargs.get('S')
            e_e = self.S / self.E
            self.e = e_e + self.e_in
            self.Ec = self.E
            # Verifica e faz a correcao plastica
            f = abs(self.S - self.beta) - self.Sy
            if f > 0:
                # Corrige a tensao e incrementa a parcela inelastica da deformacao
                self.Ec = self.Et
                r = self.S / abs(self.S)
                r2 = (self.S - self.beta) / abs(self.S - self.beta)
                C0 = (1+self.H/self.E)
                self.S = (C0*r*self.S - self.Sy - r2*self.beta) / (C0*r - r2)
                f = abs(self.S - self.beta) - self.Sy
                l = f / (self.E + self.H)
                self.e = self.e + l * r
                self.e_in = self.e_in + l * r
                self.beta = self.beta + l * self.H * r
        else:
            raise KeyError("Invalid key to function Plasticity1D.update()")

class Node:
    def __init__(self,i,x,y):
        self.i = i
        self.x = x
        self.y = y
        self.x0 = x
        self.y0 = y

class Material:
    def __init__(self,E,A,M):
        self.E = E
        self.A = A
        self.M = M

class Element:
    Rel = np.ones((4,1))
    Kel = np.ones((4,4))
    Mel = np.diag((1,1,1,1))
    def __init__(self,material,node1,node2):
        self.material = material
        self.node1 = node1
        self.node2 = node2
        dx = node2.x - node1.x
        dy = node2.y - node1.y
        self.l0 = np.sqrt(dx*dx+dy*dy)
        self.Mel = self.Mel*self.material.M*self.material.A*self.l0/2
    def update_Rloc(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        # Forcas internas
        e = (length-self.l0)/self.l0
        if isinstance(self.material.E, Plasticity1D):
            self.material.E.update(e = e)
            S = self.material.E.S
        else:
            S = self.material.E*e
        self.Rloc = S*self.material.A
    def update_Rel(self):
        self.update_Rloc()
        sen = self.sen
        cos = self.cos
        self.Rel[0][0] = -cos
        self.Rel[1][0] = -sen
        self.Rel[2][0] = cos
        self.Rel[3][0] = sen
        self.Rel = self.Rel*self.Rloc
    def update_Kloc(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        e = (length-self.l0)/self.l0
        if isinstance(self.material.E, Plasticity1D):
            self.material.E.update(e = e)
            S = self.material.E.S
            E = self.material.E.Ec
            #E = self.material.E.E
        else:
            E = self.material.E
            S = self.material.E*e
        self.Kloc = E*self.material.A/self.l0
        self.Cnl = S*self.material.A/self.length
        
    def update_Kel(self):
        ##self.update_Kloc()
        
        # Parcela Linear
        sen = self.sen
        cos = self.cos
        self.Kel[0][0] = cos*cos
        self.Kel[0][1] = cos*sen
        self.Kel[0][2] = -cos*cos
        self.Kel[0][3] = -cos*sen
        self.Kel[1][1] = sen*sen
        self.Kel[1][2] = -cos*sen
        self.Kel[1][3] = -sen*sen
        self.Kel[2][2] = cos*cos
        self.Kel[2][3] = cos*sen
        self.Kel[3][3] = sen*sen
        ##self.Kel = self.Kel*1.0
        self.Kel = self.Kel*self.Kloc
        ###print(self.Kel)
        
        # Parcela NaoLinear
        Cnl = self.Cnl
        self.Kel[0][0] = self.Kel[0][0] + Cnl*sen*sen
        self.Kel[0][1] = self.Kel[0][1] - Cnl*cos*sen
        self.Kel[0][2] = self.Kel[0][2] - Cnl*sen*sen
        self.Kel[0][3] = self.Kel[0][3] + Cnl*cos*sen
        self.Kel[1][1] = self.Kel[1][1] + Cnl*cos*cos
        self.Kel[1][2] = self.Kel[1][2] + Cnl*cos*sen
        self.Kel[1][3] = self.Kel[1][3] - Cnl*cos*cos
        self.Kel[2][2] = self.Kel[2][2] + Cnl*sen*sen
        self.Kel[2][3] = self.Kel[2][3] - Cnl*cos*sen
        self.Kel[3][3] = self.Kel[3][3] + Cnl*cos*cos
        #print(self.Kel)
        
        # Simetria
        self.Kel[1][0] = self.Kel[0][1]
        self.Kel[2][0] = self.Kel[0][2]
        self.Kel[2][1] = self.Kel[1][2]
        self.Kel[3][0] = self.Kel[0][3]
        self.Kel[3][1] = self.Kel[1][3]
        self.Kel[3][2] = self.Kel[2][3]
        ###print(self.Kel)

###A0 = 1
###E0 = 1.0e6
###m0 = 0
###A = 1.0e-3
###E = 2.08e11
###m = 0
###Et = 2.08e9
###Sy = 2.0e9
###mat0 = Material(E0,A0,m0)
###mat1 = Material(Plasticity1D(E,Et,Sy),A,m)
###mat2 = Material(E,A,m)
###L = 5
###B = L * np.cos(np.pi/12)
###H = L * np.sin(np.pi/12)
###node1 = Node(0, 0, 0)
###node2 = Node(1, L, 0)
###node3 = Node(2, L*(1+np.cos(np.pi/12)), L*np.sin(np.pi/12))
###node4 = Node(3, L*(1+np.cos(np.pi/12)), 1+L*np.sin(np.pi/12))
###nodes = [node1,node2,node3,node4]
###el1 = Element(mat2,node1,node2)
###el2 = Element(mat1,node2,node3)
####el2 = Element(mat2,node2,node3)
###el3 = Element(mat0,node3,node4)
###els = [el1,el2,el3]
###elimDOFs = [0,1,3,4,6,7]

A = 1.0e-3
E1 = 2.08e11
m = 0
E2 = 2.08e9
##Sy = 2.0e9
Sy = 1.98e9
mat1 = Material(E1,A,m)
mat2 = Material(Plasticity1D(E1,E2,Sy),A,m)
L = 1
node1 = Node(0,   0, 0)
node2 = Node(1,   L, 0)
node3 = Node(2, 2*L, 0)
nodes = [node1,node2,node3]
el1 = Element(mat1,node1,node2)
el2 = Element(mat2,node2,node3)
els = [el1,el2]
elimDOFs = [0,1,3,5]

F = [0,4.0e6]

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

##print("- Coordenadas iniciais dos nos:")
##for inode, node in enumerate(nodes):
##    print("No {0}: ({1}; {2})".format(inode,node.x,node.y))
##print()

Dt = 0.05
tot_t = 100.0
nsteps = int(tot_t / Dt)

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
T = 30.0
w = 2.0*np.pi/T
t = np.linspace(0.0, tot_t, nsteps+1)
FF = np.vstack([Fi*np.sin(w*t) for Fi in F])

LL = np.array([0])
ee = np.array([0])
SS = np.array([0])

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
    nsteps = 610
    progress_bar = lambda *args: None
else:
    dprint = lambda *args: None

step = 0
iters = []
STBK = True
while step < nsteps:

    step = step + 1
    
    ia = step - 1
    U[step] = U[ia]
    
    niter = 0
    Rnorm = np.inf
    
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
            el.update_Kloc()
            E = el.material.E
            if isinstance(E, Plasticity1D) and STBK:
            
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
                
                if niter == 1:
                    ##LL = np.append(LL,(el2.length-el2.l0)/el2.l0)
                    ee = np.append(ee,E.e)
                    SS = np.append(SS,E.S)
        
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
    
    #N = int(1.0*nsteps)
    N = int(1.0*nsteps-510)
    
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
