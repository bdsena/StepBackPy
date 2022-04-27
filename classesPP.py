import numpy as np

class Node:
    def __init__(self,i,x,y,a):
        self.i = i
        self.x = x
        self.y = y
        self.a = a
        self.x0 = x
        self.y0 = y
        self.a0 = a
        self.ainc = 0.

class MaterialLin:
    def __init__(self,E,A,m,I):
        self.E = E
        self.A = A
        self.m = m
        self.I = I

class MaterialPlast:
    def __init__(self,E,A,m,I,Mc):
        self.E = E
        self.A = A
        self.m = m
        self.I = I
        self.Mc = Mc

class Beam:

    Rel = np.ones((6,1))
    Kel = np.ones((6,6))
    
    def __init__(self,material,node1,node2):
        self.node1 = node1
        self.node2 = node2
        
        dx = node2.x - node1.x
        dy = node2.y - node1.y
        L = np.sqrt(dx*dx+dy*dy)
        
        m = material.m
        E = material.E
        A = material.A
        I = material.I
        self.m = m
        self.E = E
        self.A = A
        self.I = I
        #I = A**2/(4*np.pi)
        #I = A**2/12
        #print("E = {}".format(E))
        #print("A = {}".format(A))
        #print("I = {}".format(I))
        
        self.EA = E*A
        self.EI = E*I
        self.I = I
        
        lumped_mass = True
        if lumped_mass:
            CM = m*A*L/2
            CI = CM*L**2/12
            self.Mel = np.diag((CM,CM,CI,CM,CM,CI))
        else:
            C0 = m*A*L/420
            M11 =  C0*140
            M14 =  C0*70
            M22 =  C0*156
            M23 =  C0*22*L
            M25 =  C0*54
            M26 = -C0*13*L
            M33 =  C0*4*L**2
            M35 = -M26
            M36 = -C0*3*L**2
            M44 =  M11
            M55 =  C0*156#150
            M56 = -M23
            M66 =  M33
            self.Mel = C0*np.array([
                [M11,   0,   0, M14,   0,   0],
                [  0, M22, M23,   0, M25, M26],
                [  0, M23, M33,   0, M35, M36],
                [M14,   0,   0, M44,   0,   0],
                [  0, M25, M35,   0, M55, M56],
                [  0, M26, M36,   0, M56, M66]
            ])
        
        self.L0 = L
        cos = dx/L
        sen = dy/L
        self.theta0 = np.arctan2(sen, cos)
        self.a1 = self.node1.a0 - self.theta0
        self.a2 = self.node2.a0 - self.theta0
        
        a1cos = np.cos(self.node1.a0)
        a1sen = np.sin(self.node1.a0)
        self.a1_T = np.array([
            [ a1cos, a1sen],
            [-a1sen, a1cos],
        ])
        a2cos = np.cos(self.node2.a0)
        a2sen = np.sin(self.node2.a0)
        self.a2_T = np.array([
            [ a2cos, a2sen],
            [-a2sen, a2cos],
        ])
    
    def update_length(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        a1inc_cos = np.cos(self.node1.ainc)
        a1inc_sen = np.sin(self.node1.ainc)
        a1inc_T = np.array([
            [ a1inc_cos, a1inc_sen],
            [-a1inc_sen, a1inc_cos],
        ])
        a2inc_cos = np.cos(self.node2.ainc)
        a2inc_sen = np.sin(self.node2.ainc)
        a2inc_T = np.array([
            [ a2inc_cos, a2inc_sen],
            [-a2inc_sen, a2inc_cos],
        ])
        
        #self.a1_T = np.matmul(a1inc_T, self.a1_T)
        #self.a2_T = np.matmul(a2inc_T, self.a2_T)
        self.a1_T = np.matmul(self.a1_T, a1inc_T)
        self.a2_T = np.matmul(self.a2_T, a2inc_T)
        
        #print('self.a1_T')
        #print(self.a1_T)
        #print('self.a2_T')
        #print(self.a2_T)
        
        R = np.array([
            [ self.cos, self.sen],
            [-self.sen, self.cos],
        ])
        
        ##self.a1 = np.arcsin(.5*( - np.dot(a1inc_T[1],R[0]) + np.dot(a1inc_T[0],R[1]) ))
        ##self.a2 = np.arcsin(.5*( - np.dot(a2inc_T[1],R[0]) + np.dot(a2inc_T[0],R[1]) ))
        
        theta = np.arctan2(self.sen, self.cos)
        #print("theta = ",theta)
        #print("sen = ",self.sen)
        #print("cos = ",self.cos)
        #print("dl = ",self.length-self.L0)
        beta = np.arctan2(self.sen, self.cos) - self.theta0
        self.a1 = self.node1.a - beta
        self.a2 = self.node2.a - beta
        #print('length = ',length)
        ####print('beta = ',beta)
        #print('a1 = ',self.a1)
        #print('a2 = ',self.a2)
    
    def update_Rloc(self):
        # Forcas internas
        e = (self.length-self.L0)/self.L0
        N = self.EA*e
        M1 = 2.*self.EI/self.L0*(2.*self.a1+self.a2)
        M2 = 2.*self.EI/self.L0*(self.a1+2.*self.a2)
        V = (M1+M2)/self.L0#length
        self.Rloc = np.array([
            [-N],
            [ V],
            [M1],
            [ N],
            [-V],
            [M2],
        ])
        #print('Rloc')
        #print(self.Rloc)
    
    def update_Rel(self):
        sen = self.sen
        cos = self.cos
        T = np.array([
            [ cos, sen, 0,   0,   0, 0],
            [-sen, cos, 0,   0,   0, 0],
            [   0,   0, 1,   0,   0, 0],
            [   0,   0, 0, cos, sen, 0],
            [   0,   0, 0,-sen, cos, 0],
            [   0,   0, 0,   0,   0, 1],
        ])
        #print('T')
        #print(T)
        #print('T.T')
        #print(T.T)
        self.Rel = np.matmul(T.T, self.Rloc)
        #print('Rel')
        #print(self.Rel)
    
    def update_Kloc(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        EA = self.EA
        EI = self.EI
        L0 = self.L0
        I = self.I
        A = self.A
        L = self.length
        
        K11 = EA/L0
        K22 = 12*EI/L0**3
        K23 = 6*EI/L0**2
        K33 = 4*EI/L0
        KK33 = K33/2
        
        Kloc = np.array([
            [ K11,   0,   0,-K11,   0,   0],
            [   0, K22, K23,   0,-K22, K23],
            [   0, K23, K33,   0,-K23,KK33],
            [-K11,   0,   0, K11,   0,   0],
            [   0,-K22,-K23,   0, K22,-K23],
            [   0, K23,KK33,   0,-K23, K33],
        ])
        
        M1 = 2.*EI/L0*(2.*self.a1+self.a2)
        M2 = 2.*EI/L0*(self.a1+2.*self.a2)
        V = (M1+M2)/L0
        
        e = (self.length-L0)/L0
        F1 = -EA*e
        F2 = V
        F3 = M1#(M1+M2)/2
        ##print(self.node1.i,self.node2.i,F1,F2,F3)
        anflex = True
        if anflex:
            IAL = I/A/L0
            C1 = (1.2+12*IAL/L0)/L0
            C2 = (0.1+6*IAL/L0)
            C3 = (L0/7.5+4*IAL)
            C4 = (L0/30.-2*IAL)
            G11 = -F1/L0
            G12 = -F2/L0
            G13 = -F3/L0
            G14 = -G11
            G15 = -G12
            G16 =  F3/L0-F2
            G22 = -C1*F1
            G23 = -C2*F1
            G24 =  G15
            G25 = -G22
            G26 =  G23
            G33 = -C3*F1
            G34 = -G13
            G35 = -G23
            G36 =  C4*F1
            G44 =  G11
            G45 =  G12
            G46 = -G16
            G55 =  G22
            G56 = -G23
            G66 =  G33
            Gloc = np.array([
                [G11, G12, G13, G14, G15, G16],
                [G12, G22, G23, G24, G25, G26],
                [G13, G23, G33, G34, G35, G36],
                [G14, G24, G34, G44, G45, G46],
                [G15, G25, G35, G45, G55, G56],
                [G16, G26, G36, G46, G56, G66],
            ])
        else:
            ###C1 = F1/L0
            ###C3 = F3/L0
            ###IAL = I/A/L0
            ###C4 = F1*(1.2+12*IAL/L0)/L0
            ###C6 = F1*(0.1+6*IAL/L0)
            ###C8 = F1*(L0/7.5+4*IAL)
            ###C9 = F1*(L0/30.-2*IAL)
            C1 = 0
            C3 = 0
            C4 = F1*1.2/L0
            C6 = F1*0.1
            C8 = F1*L0/7.5
            C9 = F1*L0/30.
            Gloc = np.array([
                [  -C1,   0, -C3,    C1,   0, C3-F2],
                [    0, -C4, -C6,     0,  C4,   -C6],
                [  -C3, -C6, -C8,    C3,  C6,    C9],
                [   C1,   0,  C3,   -C1,   0, F2-C3],
                [    0,  C4,  C6,     0, -C4,    C6],
                [C3-F2, -C6,  C9, F2-C3,  C6,   -C8],
            ])
        
        self.Kloc = Kloc + Gloc
        ##print(self.Kloc)
        
    def update_Kel(self):
        sen = self.sen
        cos = self.cos
        T = np.array([
            [ cos, sen, 0,   0,   0, 0],
            [-sen, cos, 0,   0,   0, 0],
            [   0,   0, 1,   0,   0, 0],
            [   0,   0, 0, cos, sen, 0],
            [   0,   0, 0,-sen, cos, 0],
            [   0,   0, 0,   0,   0, 1],
        ])
        self.Kel = np.matmul(T.T, np.matmul(self.Kloc, T))

class HysBeam(Beam):

    f = 1
    
    def __init__(self,material,node1,node2):
        super().__init__(material,node1,node2)
    
        #Dados do material
        EI = self.EI
        EIt = EI/100
        #self.EI = EI    # Modulo de elasticidade
        self.EIt = EIt  # Modulo de elasticidade tangente
        self.Mc = material.Mc  # Momento critico
        self.H = EI * EIt / (EI - EIt) # Modulo de encruamento
    
        #Variaveis de estado
        self.k_in = 0  # Parcela inelastica das deformacoes
        self.beta = 0  # Variavel de controle da plasticidade
        self.EIc = EI  # Modulo de elasticidade atual
        self.M = 0     # Momento atual

    # Atualiza o estado de tensao
    def update_EI(self, k):
    
        # Calcula como se fosse regime elastico
        k_e = k - self.k_in
        self.M = self.EI * k_e
        self.EIc = self.EI
    
        # Verifica e faz a correcao plastica
        f = abs(self.M - self.beta) - self.Mc
        if (f > 0):
            # Corrige a tensao e incrementa a parcela inelastica da deformacao
            self.EIc = self.EIt
            r = self.M / abs(self.M)
            l = f / (self.EI + self.H)
            self.M = self.M - self.EI * l * r
            self.k_in = self.k_in + l * r
            self.beta = self.beta + l * self.H * r
    
    def update_length(self):
        super().update_length()
        k = (self.a1 + self.a2)/2
        self.update_EI(k)
    
    def update_Kloc(self):
        #EI_bak = self.EI
        #self.EI = self.EIc
        super().update_Kloc()
        #self.EI = EI_bak
