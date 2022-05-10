from classesPP import Node, MaterialLin, MaterialPlast, Beam, HysBeam

A = 5.0669962511586369e-004
I = 2.0431712321000111e-008
#A = 6.4516e-004
#I = 3.4686e-008
E1 = 2.06844e8
m = 0
##E2 = 2.08e9
##E2 = 3.08e10
##Sy = 2.0e9
##Sy = 1.98e9
mat1 = MaterialLin(E1,A,m,I)
mat2 = MaterialPlast(E1,A,m,I,2.,E1*I/3.)
L = 2.54
npts = 11
nodes = []
for i in range(npts):
    nodes.append(Node(i, L*i/(npts-1), 0, 0))
els = []
for i in range(npts-1):
    #els.append(Beam(i,mat1,nodes[i],nodes[i+1]))
    els.append(HysBeam(i,mat2,nodes[i],nodes[i+1]))

Upresc = {
    0: 0.0, # no' 0 x
    1: 0.0, # no' 0 y
    2: 0.0, # no' 0 a
    #30: 0.0, # no' 11 x
    31: 1.0, # no' 11 y
}
Mc = 0.5*np.pi*els[0].EI/L
Fc = els[0].EI/L**2
Fpresc = {
    #31: Fc*2, # no' 11 y
    #32: Mc, # no' 11 a
}
print('Upresc')
print(Upresc)
print('Fpresc')
print(Fpresc)

T = 20.0

STBK = False

Dt = 1.
#tot_t = 3.55
#tot_t = 6.3
tot_t = 40.
