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
            dprint(self.S,self.beta,self.Sy)
            f = abs(self.S - self.beta) - self.Sy
            if f > 0:
                # Corrige a tensao e incrementa a parcela inelastica da deformacao
                self.Ec = self.Et
                r = self.S / abs(self.S)
                l = f / (self.E + self.H)
                self.S = self.S - self.E * l * r
                self.e_in = self.e_in + l * r
                self.beta = self.beta + l * self.H * r
                dprint("S = {}".format(self.S))
                dprint("e_e = {}".format(e_e))
                dprint("e = {}".format(self.e))
                dprint("f = {}".format(f))
                dprint("r = {}".format(r))
                dprint("l = {}".format(l))
                dprint("e_in = {}".format(self.e_in))
                dprint("beta = {}".format(self.beta))
        elif 'S' in kwargs:
            # Calcula como se fosse regime elastico
            self.S = kwargs.get('S')
            e_e = self.S / self.E
            self.e = e_e + self.e_in
            self.Ec = self.E
            # Verifica e faz a correcao plastica
            dprint(self.S,self.beta,self.Sy)
            f = abs(self.S - self.beta) - self.Sy
            if f > 0:
                # Corrige a tensao e incrementa a parcela inelastica da deformacao
                self.Ec = self.Et
                r = self.S / abs(self.S)
                r2 = (self.S - self.beta) / abs(self.S - self.beta)
                C0 = (1+self.H/self.E)
                dprint(self.S,C0,r,r2,self.beta)
                self.S = (C0*r*self.S - self.Sy - r2*self.beta) / (C0*r - r2)
                f = abs(self.S - self.beta) - self.Sy
                l = f / (self.E + self.H)
                self.e = self.e + l * r
                self.e_in = self.e_in + l * r
                self.beta = self.beta + l * self.H * r
                dprint("S = {}".format(self.S))
                dprint("e_e = {}".format(e_e))
                dprint("e = {}".format(self.e))
                dprint("f = {}".format(f))
                dprint("r = {}".format(r))
                dprint("r2 = {}".format(r2))
                dprint("l = {}".format(l))
                dprint("e = {}".format(self.e))
                dprint("e_in = {}".format(self.e_in))
                dprint("beta = {}".format(self.beta))
        else:
            raise KeyError("Invalid key to function Plasticity1D.update()")

E1 = 2.08e11
E2 = 3.08e10
Sy = 1.98e9
EE = Plasticity1D(E1,E2,Sy)

Dt = 0.05
tot_t = 3.55#6.3#20.0#100.0
nsteps = int(tot_t / Dt)

t = np.linspace(0.0, tot_t, nsteps+1)
N = int(1.0*nsteps)

def dprint(*args):
    print(*args)

DEBUG = False
if DEBUG:
    N = 11
else:
    dprint = lambda *args: None

data1 = np.load('out.npz')
data2 = np.load('outSB.npz')
SS1 = data1['SS'][1]
ee1 = data1['ee'][1]
SS2 = data2['SS'][1]
ee2 = data2['ee'][1]

step = 0
SS3 = np.copy(SS2)
ee3 = np.array([0])
while step < N:
    step = step + 1
    print("## step {}".format(step))
    print("S = ",SS2[step])
    print("e_in = ", EE.e_in)
    print("update")
    EE.update(S = SS2[step])
    print("e_in = ", EE.e_in)
    ee3 = np.append(ee3,EE.e)

exit()
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,5))
#plt.subplot(311)
#plt.plot(t[:N],SS1[:N],'-')
#plt.plot(t[:N],SS2[:N],'-')
#plt.xlabel("t")
#plt.ylabel("S")
#plt.subplot(312)
#plt.plot(t[:N],ee1[:N],'-')
#plt.plot(t[:N],ee2[:N],'-')
#plt.xlabel("t")
#plt.ylabel("e")
#plt.subplot(313)
plt.plot(ee1[-3:],SS1[-3:],'-o')
plt.plot(ee2[-3:],SS2[-3:],'-o')
plt.plot(ee3[-3:],SS3[-3:],'-o')
plt.xlabel("e")
plt.ylabel("S")
plt.show()
