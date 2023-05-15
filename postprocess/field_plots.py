#-------------------------------------- Libraries -----------------------------------------#
import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from pylab import rcParams
mpl.rcParams.update(mpl.rcParamsDefault)
#------------------------------------------------------------------------------------------#

#----------------------------------- Import other files -----------------------------------#
import plot as pl
import sys
sys.path.append("..")
import para
sys.path.append("../lib")
import grid
from derivative import *
#------------------------------------------------------------------------------------------#

mpl.style.use('classic')

if pl.type == "field":
    i = pl.t_start
    j = 0
    while i <= para.tfinal:
        with h5py.File(pl.path + pl.output_folder + "/fields/2D_%d.00.h5" %(i), "r") as f:
            # List all groups
            # print("Keys: %s" % f.keys())
            
            data_T = list(f.keys())[0]
            data_rho = list(f.keys())[2]
            data_ux = list(f.keys())[3]
            data_uz = list(f.keys())[4]
            data_z_p = list(f.keys())[1]

            # Get the data
            T = np.array(f[data_T])
            rho = np.array(f[data_rho])
            ux = np.array(f[data_ux])
            uz = np.array(f[data_uz])
            param = np.array(f[data_z_p])

        fig = plt.figure(figsize=(8,4), tight_layout = True)
        
        # For temperature
        ax1 = fig.add_subplot(121)
        c1 = ax1.pcolor(grid.X_mesh, grid.Z_mesh, T, cmap=cm.coolwarm)
        ax1.set_aspect(aspect = 1)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb1 = fig.colorbar(c1, cax=cax)
        cb1.ax.tick_params(labelsize=10)
        # cb1.set_ticks([1,1+0.5*para.epsilon,1+para.epsilon])
        ax1.quiver(grid.X_mesh[::pl.Nq,::pl.Nq],grid.Z_mesh[::pl.Nq,::pl.Nq],ux[::pl.Nq,::pl.Nq],uz[::pl.Nq,::pl.Nq],color = 'k',label='$\vec{u}$', scale=0.15, scale_units='inches')
        ax1.set_xlim(0,para.A)
        ax1.set_ylim(0,1)
        ax1.tick_params(axis='both', which='major', labelsize=10)
        ax1.set_xlabel(r'$X$', fontsize=10)
        ax1.set_ylabel(r'$Z$', fontsize=10)
        ax1.set_title(r'$T$', fontsize=12)
        
        # For density
        ax2 = fig.add_subplot(122)
        c2 = ax2.pcolor(grid.X_mesh, grid.Z_mesh, rho, cmap=cm.coolwarm)
        ax2.set_aspect(aspect = 1)
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb2 = fig.colorbar(c2, cax=cax)
        cb2.ax.tick_params(labelsize=10)
        # cb2.set_ticks([1,np.round(1+0.5*grid.Chi,1),grid.Chi])
        ax2.quiver(grid.X_mesh[::pl.Nq,::pl.Nq],grid.Z_mesh[::pl.Nq,::pl.Nq],ux[::pl.Nq,::pl.Nq],uz[::pl.Nq,::pl.Nq],color = 'k',label='$\vec{u}$', scale=0.15, scale_units='inches')
        ax2.set_xlim(0,para.A)
        ax2.set_ylim(0,1)
        ax2.tick_params(axis='both', which='major', labelsize=10)
        ax2.set_xlabel(r'$X$', fontsize=10)
        ax2.set_ylabel(r'$Z$', fontsize=10)
        ax2.set_title(r'$\rho$', fontsize=12)
        
        fig.suptitle(r'$\alpha = %.1f$,    ' %(param[0]) + r'$m = %.1f$,    ' %(param[1]) + r'$\epsilon = %d$,    ' %(para.epsilon) + r'$Pr = %d$,    ' %(param[2]) + r'$Ra = %d$,    ' %(param[3]) + r'$t = %d$' %(grid.t_f_step[j]), fontsize=15)  
        # plt.rc('text', usetex=True) 
        plt.savefig(pl.path + pl.output_folder + '/plots/field_plot/%d.00.png' %(i))
        # plt.show()

        i = i+para.t_f
        j = j+1
        pass
    pass

