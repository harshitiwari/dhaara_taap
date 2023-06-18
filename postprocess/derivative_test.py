import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from pylab import rcParams

import sys
sys.path.append("..")
import para

sys.path.append("../lib")
import grid
from derivative import *

# Some analytical function
A = np.sin(grid.X_mesh)*np.cos(grid.Z_mesh)

dA_x_ana = np.cos(grid.X_mesh)*np.cos(grid.Z_mesh)      # Analytical derivative w.r.t x
dA_z_ana = -np.sin(grid.X_mesh)*np.sin(grid.Z_mesh)     # Analytical derivative w.r.t z

# Errors
del_x = dA_x_ana - dfx_c(A)
del_z = dA_z_ana - dfz_c(A)

# For polar gird
grad_p = dfz_c(grid.T_0*grid.rho_0)
g_term = -grid.rho_0*grid.C3/grid.C1

# Plotting errors
fig = plt.figure()

if para.coordinates == 'cartesian':

    ax1 = fig.add_subplot(121,aspect = 'equal')
    c1 = ax1.pcolormesh(grid.X_mesh,grid.Z_mesh,del_x, shading='auto')
    ax1.set_aspect(aspect = 1)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(c1, cax=cax1)
    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    ax1.set_xlabel('$X$')
    ax1.set_ylabel('$Z$')
    ax1.set_title(r'$\Delta_x$')

    ax2 = fig.add_subplot(122,aspect = 'equal')
    c2 = ax2.pcolormesh(grid.X_mesh,grid.Z_mesh,del_z, shading='auto')
    ax2.set_aspect(aspect = 1)
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(c2, cax=cax2)
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_xlabel('$X$')
    ax2.set_ylabel('$Z$')
    ax2.set_title(r'$\Delta_z$')

    plt.show()

elif para.coordinates == 'polar':

    ax1 = fig.add_subplot()
    c1 = ax1.pcolormesh(grid.X_mesh,grid.Z_mesh,grad_p - g_term, shading='auto')
    ax1.set_aspect(aspect = 1)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(c1, cax=cax1)
    ax1.set_xlim(0,2*np.pi)
    ax1.set_ylim(para.r_b,1)
    ax1.set_xlabel('$X$')
    ax1.set_ylabel('$Z$')
    ax1.set_title('Error')
    plt.show()