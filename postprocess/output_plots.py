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

data = np.genfromtxt(pl.path + pl.output_folder + pl.output_file)

# First two lines contain parameters and last line contain time of simulation
t = data[:,0]
Ke = data[:,1]
Ie = data[:,2]
V_rms = data[:,3]
F_c = data[:,4]
F_r = data[:,5]
F_k = data[:,6]
cri = data[:,7]
p_t = data[:,8]

if pl.type == 'energy':
    plt.figure(figsize=(5,3))
    plt.loglog(t, Ke, label = r'$K_e = \frac{1}{2}\rho V^2$', color = 'r')
    plt.loglog(t, Ie, label = r'$I_e = \rho T$', color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$E$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim(0,1)
    plt.xticks([1,1e1,1e2],fontsize=10)
    plt.yticks([1e-6,1e-3,1e1],fontsize=10)
    plt.legend(loc=3, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/energy.png', dpi =200)

elif pl.type == 'Nusselt':
    plt.figure(figsize=(5,3))
    plt.plot(t, (F_r)/(-grid.C4*para.epsilon), label = r'$\mathrm{Nu}$', color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$\mathrm{Nu}$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim(1,1)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc=4, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/Nusselt.png', dpi =200)

elif pl.type == 'criteria':
    plt.figure(figsize=(5,3))
    plt.plot(t, np.absolute(cri), label = r'$\langle \left| \frac{\partial T}{\partial z} \right| - \frac{g}{C_p} \rangle$', color = 'k')
    plt.xlabel(r'$t$', fontsize=13)           
    plt.xlim(0, para.tfinal)      
    # plt.ylim(4e-2,4e-4+4e-2)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc=4, prop={'size': 18}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/criteria.png', dpi =200)

elif pl.type == 'ratio':
    plt.figure(figsize=(5,3))
    plt.plot(t[::100], (-p_t[::100]/grid.C3), label = r'$-\langle \frac{1}{\rho} \frac{\partial p}{\partial z} \rangle / g $', color = 'k')
    # plt.plot(t, -g_t, label = r'$\langle \rho g \rangle$', color = 'r')
    plt.axhline(y = 1, color = 'k', linestyle = '--', dashes=(3, 2), linewidth=1)
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim(-1.5e-4+1,1e-4+1)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc="center right", prop={'size': 18}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/ratio.png', dpi =200)

elif pl.type == 'V_rms':
    plt.figure(figsize=(5,3))
    plt.plot(t, V_rms, color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$\langle V_{rms} \rangle$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    plt.plot([], [], ' ', label=r'$\mathrm{Re} = %.2f$' %(V_rms[-1]/grid.C2))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim(0,0.8e-1)
    plt.xticks(fontsize=10)
    plt.yticks([0,0.2e-1,0.4e-1,0.6e-1,0.8e-1],fontsize=10)
    plt.legend(loc=4, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/V_rms.png', dpi =200)