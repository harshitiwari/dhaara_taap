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

data =np.genfromtxt(pl.path + pl.output_folder + pl.output_file)

# First two lines contain parameters and last line contain time of simulation
t = data[:,0]
Ke = data[:,1]
Ie = data[:,2]
V_rms = data[:,3]
Nu_b = data[:,4]
Nu_t = data[:,5]
cri = data[:,6]
p_t = data[:,7]
g_t = data[:,8]

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
    plt.savefig('plots/energy.png')

elif pl.type == 'Nusselt':
    plt.figure(figsize=(5,3))
    plt.plot(t, Nu_b, label = r'$\mathrm{Nu}_b$', color = 'k')
    plt.plot(t, Nu_t, label = r'$\mathrm{Nu}_t$', color = 'r')
    plt.xlabel(r'$t$', fontsize=13)        
    plt.ylabel(r'$\mathrm{Nu}$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim(0,1)
    plt.xticks(fontsize=10)
    plt.yticks([1,1.1,1.2,1.3],fontsize=10)
    plt.legend(loc=1, prop={'size': 15}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig('plots/Nusselt.png')

elif pl.type == 'criteria':
    plt.figure(figsize=(5,3))
    plt.plot(t, np.absolute(cri), label = r'$\langle \left| \frac{\partial T}{\partial z} \right| - \frac{g}{C_p} \rangle$', color = 'k')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$Nu$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim(0,1)
    plt.xticks(fontsize=10)
    plt.yticks([0.405,0.41,0.415,0.42], fontsize=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc=4, prop={'size': 18}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig('plots/criteria.png')

elif pl.type == 'ratio':
    plt.figure(figsize=(5,3))
    plt.plot(t, (p_t+g_t), label = r'$\langle \frac{\partial p}{\partial z} \rangle + \langle \rho g \rangle$', color = 'k')
    # plt.plot(t, -g_t, label = r'$\langle \rho g \rangle$', color = 'r')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13)     
    plt.xlim(0, para.tfinal)      
    # plt.ylim()
    plt.xticks(fontsize=10)
    plt.yticks([-0.5e-4,0,1e-4,2e-4],fontsize=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc=4, prop={'size': 18}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig('plots/ratio.png')

