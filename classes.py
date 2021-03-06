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

    # Atualiza o estado de tensao
    def update(self, e):
    
        # Calcula como se fosse regime elastico
        e_e = e - self.e_in
        self.S = self.E * e_e
        self.Ec = self.E
    
        # Verifica e faz a correcao plastica
        f = abs(self.S - self.beta) - self.Sy
        if (f > 0):
            # Corrige a tensao e incrementa a parcela inelastica da deformacao
            self.Ec = self.Et
            r = self.S / abs(self.S)
            l = f / (self.E + self.H)
            self.S = self.S - self.E * l * r
            self.e_in = self.e_in + l * r
            self.beta = self.beta + l * self.H * r

class Node:
    def __init__(self,i,x,y):
        self.i = i
        self.x = x
        self.y = y
        self.x0 = x
        self.y0 = y

class MaterialLin:
    def __init__(self,E,A,M):
        self.E = E
        self.A = A
        self.M = M

class MaterialPlast:
    def __init__(self,E1,E2,Sy,A,M):
        self.E = Plasticity1D(E1,E2,Sy)
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
    
    def update_length(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
    
    def update_Rloc(self):
        # Forcas internas
        e = (self.length-self.l0)/self.l0
        S = self.material.E*e
        self.Rloc = S*self.material.A
    
    def update_Rel(self):
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
        E = self.material.E
        S = self.material.E*e
        self.Kloc = E*self.material.A/self.l0
        self.Cnl = S*self.material.A/self.length
        
    def update_Kel(self):
        
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
        
        # Simetria
        self.Kel[1][0] = self.Kel[0][1]
        self.Kel[2][0] = self.Kel[0][2]
        self.Kel[2][1] = self.Kel[1][2]
        self.Kel[3][0] = self.Kel[0][3]
        self.Kel[3][1] = self.Kel[1][3]
        self.Kel[3][2] = self.Kel[2][3]

class HysElement(Element):

    f = 1
    
    def update_Rloc(self):
        # Forcas internas
        e = (self.length-self.l0)/self.l0
        E = self.material.E
        E.update(e)
        self.f = E.Ec/E.E
        S = E.S
        self.Rloc = S*self.material.A
    
    def update_Kloc(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        e = (length-self.l0)/self.l0
        S = self.material.E.S
        E = self.material.E.E*self.f
        self.Kloc = E*self.material.A/self.l0
        self.Cnl = S*self.material.A/self.length
