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
    def __init__(self,E,A,m,I,Mc,EIt):
        self.E = E
        self.A = A
        self.m = m
        self.I = I
        self.Mc = Mc
        self.EIt = EIt

class Beam:

    Rel = np.ones((6,1))
    Kel = np.ones((6,6))
    
    def __init__(self,iel,material,node1,node2):
        self.node1 = node1
        self.node2 = node2
        self.iel = iel
        
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
    
    def set_STBK(self, STBK):
        pass
    
    def backup(self):
        pass
    
    def restore(self):
        pass
    
    def NRP_check(self):
        return False

    def update_length(self):
        dx = self.node2.x - self.node1.x
        dy = self.node2.y - self.node1.y
        length = np.sqrt(dx*dx+dy*dy)
        self.cos = dx/length
        self.sen = dy/length
        self.length = length
        
        self.theta = np.arctan2(self.sen, self.cos) - self.theta0
        self.a1 = (self.node1.a - self.theta)
        self.a2 = (self.node2.a - self.theta)
        
        self.k1 = -2.*(2.*self.a1+self.a2)/self.L0
        self.k2 =  2.*(self.a1+2.*self.a2)/self.L0
    
    def update_Rloc(self):
        e = (self.length-self.L0)/self.L0
        N = self.EA*e
        M1 = -self.EI*self.k1
        M2 =  self.EI*self.k2
        V = (M1+M2)/self.L0
        self.Rloc = np.array([
            [-N],
            [ V],
            [M1],
            [ N],
            [-V],
            [M2],
        ])
    
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
        self.Rel = np.matmul(T.T, self.Rloc)
    
    def update_KL(self):
        
        EA = self.EA
        EI = self.EI
        L0 = self.L0
        
        K11 = EA/L0
        K22 = 12*EI/L0**3
        K23 = 6*EI/L0**2
        K33 = 4*EI/L0
        KK33 = K33/2
        self.KL = np.array([
            [ K11,   0,   0,-K11,   0,   0],
            [   0, K22, K23,   0,-K22, K23],
            [   0, K23, K33,   0,-K23,KK33],
            [-K11,   0,   0, K11,   0,   0],
            [   0,-K22,-K23,   0, K22,-K23],
            [   0, K23,KK33,   0,-K23, K33],
        ])
    
    def update_KNL(self):
        
        F1 = self.Rloc[0][0]
        F2 = self.Rloc[1][0]
        F3 = self.Rloc[2][0]

        I = self.I
        A = self.A
        L0 = self.L0

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
        self.KNL = np.array([
            [G11, G12, G13, G14, G15, G16],
            [G12, G22, G23, G24, G25, G26],
            [G13, G23, G33, G34, G35, G36],
            [G14, G24, G34, G44, G45, G46],
            [G15, G25, G35, G45, G55, G56],
            [G16, G26, G36, G46, G56, G66],
        ])
    
    def update_Kloc(self, step):
        self.step = step
        self.update_KL()
        self.update_KNL()
        self.Kloc = self.KL + self.KNL
        
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
    
    f = np.ones(2)
    M = np.zeros(2)
    k = np.zeros(2)
    k_in = np.zeros(2)
    beta = np.zeros(2)
    M_bak = np.zeros(2)
    k_bak = np.zeros(2)
    k_in_bak = np.zeros(2)
    beta_bak = np.zeros(2)

    def __init__(self,iel,material,node1,node2):
        super().__init__(iel,material,node1,node2)
        self.Mc = material.Mc
        self.EIt = material.EIt
        self.H = self.EI * self.EIt / (self.EI - self.EIt)
    
    def set_STBK(self, STBK):
        self.STBK = STBK

    def backup(self):
        self.k_in_bak = self.k_in
        self.beta_bak = self.beta
        self.M_bak = self.M
        self.k_bak = self.k

    def restore(self):
        self.k_in = self.k_in_bak
        self.beta = self.beta_bak
        self.M = self.M_bak
        self.k = self.k_bak
    
    def NRP_check(self):
        
        NRP_Loop = False
        
        if self.STBK:
            self.update_k()
            dM = self.M - self.M_bak
            dk = self.k - self.k_bak
            EIsec = abs(dM/dk)
            EIsec = np.where(np.isnan(EIsec),0,EIsec)
            f_new = EIsec/self.EI
            err = (f_new - self.f) / self.f
            if np.linalg.norm(err) > 1e-2:
                self.f = f_new
                NRP_Loop = True
        
        return NRP_Loop
    
    def update_k(self):
        # Calcula como se fosse regime elastico
        M = self.M
        k_e = M / self.EI
        self.k = k_e + self.k_in_bak
        # Verifica e faz a correcao plastica
        f = abs(M - self.beta_bak) - self.Mc
        r = M / abs(M)
        r2 = (M - self.beta_bak) / abs(M - self.beta_bak)
        C0 = (1+self.H/self.EI)
        M = (C0*r*M - self.Mc - r2*self.beta_bak) / (C0*r - r2)
        f2 = abs(M - self.beta_bak) - self.Mc
        l = f2 / (self.EI + self.H)
        lr = np.where(f > 0, l*r, 0.)
        self.k = self.k + lr
        self.k_in = self.k_in_bak + lr
        self.beta = self.beta_bak + lr * self.H
                
    def update_M(self):
        if self.STBK:
            dk = self.k - self.k_bak
            EIc = self.f * self.EI
            self.M = self.M_bak + EIc*dk
        else:
            # Calcula como se fosse regime elastico
            k_e = self.k - self.k_in_bak
            self.M = self.EI * k_e
            # Verifica e faz a correcao plastica
            f = abs(self.M - self.beta_bak) - self.Mc
            r = self.M / abs(self.M)
            l = f / (self.EI + self.H)
            lr = np.where(f > 0, l*r, 0.)
            self.M = self.M - self.EI * lr
            self.k_in = self.k_in_bak + lr
            self.beta = self.beta_bak + lr * self.H
    
    def update_length(self):
        super().update_length()
        self.k = np.array([self.k1, self.k2])

    def update_Rloc(self):

        super().update_Rloc()
        
        self.update_M()
        M1 = -self.M[0]
        M2 =  self.M[1]
        V = (M1+M2)/self.L0
        self.Rloc[1] = V
        self.Rloc[2] = M1
        self.Rloc[4] = -V
        self.Rloc[5] = M2

        ##self.k = (self.a2-self.a1)/self.L0
        ###self.k = (self.k1-self.k2)/2
        ##self.update_M()
        ##if abs(self.k) > 1e-5:
        ##    self.Rloc[2] = self.M*self.k1/self.k
        ##    self.Rloc[5] = self.M*self.k2/self.k
        ##else:
        ##    self.Rloc[2] = 0.0
        ##    self.Rloc[5] = 0.0
    
    def update_KL(self):
        
        #super().update_KL()

        EA = self.EA
        EI = self.EI
        L0 = self.L0
        f1 = self.f[0]
        f2 = self.f[1]
        
        K11 = EA/L0
        K22 = 12*EI/L0**3
        K23 = 6*EI/L0**2
        K33 = 4*EI/L0
        KK33 = K33/2
        self.KL = np.array([
            [ K11,   0*f1,   0*f1,-K11,   0*f2,   0*f2],
            [   0, K22*f1, K23*f1,   0,-K22*f2, K23*f2],
            [   0, K23*f1, K33*f1,   0,-K23*f2,KK33*f2],
            [-K11,   0*f1,   0*f1, K11,   0*f2,   0*f2],
            [   0,-K22*f1,-K23*f1,   0, K22*f2,-K23*f2],
            [   0, K23*f1,KK33*f1,   0,-K23*f2, K33*f2],
        ])

        ###EA = self.EA
        ###EI = self.EI
        ###L0 = self.L0
        ###
        ###K11 = EA/L0
        ###K22 = 12*EI/L0**3
        ###K23 = 6*EI/L0**2
        ###K33 = 4*EI/L0
        ###KK33 = K33/2
        ###M1 = self.Rloc[2][0]
        ###M2 = self.Rloc[5][0]
        ###if abs(self.a1) > 0. and abs(self.a2) > 0.:
        ###    x1 = self.L0*(2*M1-M2)/(6*self.a1*EI)
        ###    x2 = self.L0*(2*M2-M1)/(6*self.a2*EI)
        ###else:
        ###    x1 = 1.0
        ###    x2 = 1.0
        ###self.KL = np.array([
        ###    [ K11,   0,   0*x1,-K11,   0,   0*x2],
        ###    [   0, K22, K23*x1,   0,-K22, K23*x2],
        ###    [   0, K23, K33*x1,   0,-K23,KK33*x2],
        ###    [-K11,   0,   0*x1, K11,   0,   0*x2],
        ###    [   0,-K22,-K23*x1,   0, K22,-K23*x2],
        ###    [   0, K23,KK33*x1,   0,-K23, K33*x2],
        ###])
        ###disps = np.array([
        ###    [0],
        ###    [0],
        ###    [self.a1],
        ###    [self.length-self.L0],
        ###    [0],
        ###    [self.a2],
        ###])
        ###Rloc_test = np.matmul(self.KL,disps)
        ###a = 1
