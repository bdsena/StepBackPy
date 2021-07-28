from classes import Material, Node, Element, Plasticity1D

A = 1.0e-3
E1 = 2.08e11
m = 0
##E2 = 2.08e9
E2 = 3.08e10
##Sy = 2.0e9
Sy = 1.98e9
mat1 = Material(E1,A,m)
mat2 = Material(Plasticity1D(E1,E2,Sy),A,m)
L = 1
node1 = Node(0,   0, 0)
node2 = Node(1,   L, 0)
node3 = Node(2, 2*L, 0)
##node3 = Node(2, L*(1+np.cos(np.pi/12)), L*np.sin(np.pi/12))
nodes = [node1,node2,node3]
el1 = Element(mat1,node1,node2)
el2 = Element(mat2,node2,node3)
els = [el1,el2]
elimDOFs = [0,1,3,5]
##elimDOFs = [0,1,3,4]

T = 25.07
F = [0,2.5e6]

Dt = 0.05
tot_t = 100.0