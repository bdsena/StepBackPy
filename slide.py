import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

try:
    fname = sys.argv[1]
except:
    raise ValueError()

data = np.load(fname)
t = data['t']
desl = data['desl']
FF = data['FF'][1]
SS = data['SS']
ee = data['ee']
xnodes = data['xnodes']
ynodes = data['ynodes']
dofs = data['dofs']
ndof = data['ndof']

nnodes = int(len(dofs)/ndof)
nel = len(SS)

nt = len(t)
dt = t[1]-t[0]

# Create the figure and the line that we wiee manipulate
fig, ax = plt.subplots(2,2)
line1, = ax[0,0].plot(t, FF, lw=2)
line2, = ax[0,1].plot(t, desl[0])
line3, = ax[1,0].plot(xnodes[-1], ynodes[-1], '-o')
elemsel, = ax[1,0].plot(xnodes[-1][:2], ynodes[-1][:2], '-o', color='orange')
nodesel, = ax[1,0].plot([xnodes[-1][0]], [ynodes[-1][0]], 'r-o')
line4, = ax[1,1].plot(ee[0], SS[0])
line4t, = ax[1,1].plot(ee[0][-1], SS[0][-1], 'r-o')
ax[0,0].set_xlabel('Time [s]')

# adjust the main plot to make room for the sliders
plt.subplots_adjust(bottom=0.20)
# Make a horizontal slider to control the frequency.
ax_tslider = plt.axes([0.15, 0.11, 0.75, 0.01])
ax_elslider = plt.axes([0.15, 0.01, 0.75, 0.01])
ax_nodeslider = plt.axes([0.15, 0.06, 0.40, 0.01])
ax_dofslider = plt.axes([0.65, 0.06, 0.20, 0.01])
tslider = Slider(ax=ax_tslider,label='Time (s)',valmin=t[0],valmax=t[-1],valinit=t[-1],valstep=dt)
elslider = Slider(ax=ax_elslider,label='Element',valmin=0,valmax=nel-1,valinit=0,valstep=1)
nodeslider = Slider(ax=ax_nodeslider,label='Node',valmin=0,valmax=nnodes-1,valinit=0,valstep=1)
dofslider = Slider(ax=ax_dofslider,label='DOF',valmin=0,valmax=ndof-1,valinit=0,valstep=1)

def resize_x(ax, xdata):
    xmin = xdata.min()
    xmax = xdata.max()
    rang = xmax-xmin
    fac = 0.05
    ax.set_xlim(xmin=xmin-fac*rang, xmax=xmax+fac*rang)
def resize_y(ax, ydata):
    ymin = ydata.min()
    ymax = ydata.max()
    rang = ymax-ymin
    fac = 0.05
    ax.set_ylim(ymin=ymin-fac*rang, ymax=ymax+fac*rang)

def update(val):
    npt = int(tslider.val/dt)
    iel = int(elslider.val)
    inode = int(nodeslider.val)
    idof = int(dofslider.val)
    idesl = dofs[inode*ndof+idof]
    line1.set_data(t[:npt], FF[:npt])
    if np.isnan(idesl):
        line2.set_data([], [])
    else:
        ydata = desl[int(idesl)]
        line2.set_data(t[:npt], ydata[:npt])
        resize_y(ax[0,1], ydata)
    line3.set_data(xnodes[npt], ynodes[npt])
    nodesel.set_data(xnodes[npt][inode], ynodes[npt][inode])
    elemsel.set_data(xnodes[npt][iel:iel+2], ynodes[npt][iel:iel+2])
    line4.set_data(ee[iel][:npt], SS[iel][:npt])
    line4t.set_data(ee[iel][npt], SS[iel][npt])
    resize_x(ax[1,1], ee[iel])
    resize_y(ax[1,1], SS[iel])
    fig.canvas.draw_idle()

tslider.on_changed(update)
elslider.on_changed(update)
nodeslider.on_changed(update)
dofslider.on_changed(update)

plt.show()
