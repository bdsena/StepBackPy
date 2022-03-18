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
##node3 = Node(2, 2*L, 0)
ang = np.pi/12
node3 = Node(2, L*(1+np.cos(ang)), L*np.sin(ang))
nodes = [node1,node2,node3]
el1 = Element(mat1,node1,node2)
el2 = Element(mat2,node2,node3)
els = [el1,el2]
elimDOFs = [0,1,3,5]
##elimDOFs = [0,1,3,4]

T = 25.07
F = [0,2.5e6]

Dt = 0.05
#tot_t = 3.55
#tot_t = 6.3
tot_t = 100.0

#############A0 = 1
#############E0 = 1.0e6
#############m0 = 0
#############A = 1.0e-3
#############E = 2.08e11
#############m = 0
#############Et = 2.08e9
#############Sy = 2.0e9
#############mat0 = Material(E0,A0,m0)
#############mat1 = Material(Plasticity1D(E,Et,Sy),A,m)
#############mat2 = Material(E,A,m)
#############L = 5
#############B = L * np.cos(np.pi/12)
#############H = L * np.sin(np.pi/12)
#############node1 = Node(0, 0, 0)
#############node2 = Node(1, L, 0)
#############node3 = Node(2, L*(1+np.cos(np.pi/12)), L*np.sin(np.pi/12))
#############node4 = Node(3, L*(1+np.cos(np.pi/12)), 1+L*np.sin(np.pi/12))
#############nodes = [node1,node2,node3,node4]
#############el1 = Element(mat2,node1,node2)
#############el2 = Element(mat1,node2,node3)
##############el2 = Element(mat2,node2,node3)
#############el3 = Element(mat0,node3,node4)
#############els = [el1,el2,el3]
#############elimDOFs = [0,1,3,4,6,7]
#############
#############F = [0,4.0e6]
#############T = 37.77
#############
#############Dt = 0.05
#############tot_t = 100.0
#############