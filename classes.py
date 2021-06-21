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