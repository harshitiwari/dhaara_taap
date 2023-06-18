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
import output_plots as op
#------------------------------------------------------------------------------------------#

mpl.style.use('classic')

with h5py.File(pl.path + pl.output_folder + '/timeseries.h5', "r") as f:
    # List all groups
    # print("Keys: %s" % f.keys())

    data_nT_1 = list(f.keys())[0]
    data_nT_1_del = list(f.keys())[1]
    data_nT_2 = list(f.keys())[2]
    data_nT_2_del = list(f.keys())[3]
    data_nux = list(f.keys())[4]
    data_nuz = list(f.keys())[5]
    data_px = list(f.keys())[6]
    data_px_0 = list(f.keys())[7]
    data_pz = list(f.keys())[8]
    data_pz_0 = list(f.keys())[9]

    nT_1 = np.array(f[data_nT_1])
    nT_2 = np.array(f[data_nT_2])
    nT_1_del = np.array(f[data_nT_1_del])
    nT_2_del = np.array(f[data_nT_2_del])
    nux = np.array(f[data_nux])
    nuz = np.array(f[data_nuz])
    px = np.array(f[data_px])
    pz = np.array(f[data_pz])
    px_0 = np.array(f[data_px_0])
    pz_0 = np.array(f[data_pz_0])

    nT_1_i, nT_1_ii, nT_1_iii, nT_1_iv = nT_1[:,0], nT_1[:,1], nT_1[:,2], nT_1[:,3]
    nT_2_i, nT_2_ii, nT_2_iii, nT_2_iv = nT_2[:,0], nT_2[:,1], nT_2[:,2], nT_2[:,3]
    nT_1_del_i, nT_1_del_ii, nT_1_del_iii, nT_1_del_iv = nT_1_del[:,0], nT_1_del[:,1], nT_1_del[:,2], nT_1_del[:,3]
    nT_2_del_i, nT_2_del_ii, nT_2_del_iii, nT_2_del_iv = nT_2_del[:,0], nT_2_del[:,1], nT_2_del[:,2], nT_2_del[:,3]
    nux_i, nux_ii, nux_iii, nux_iv = nux[:,0], nux[:,1], nux[:,2], nux[:,3]
    nuz_i, nuz_ii, nuz_iii, nuz_iv = nuz[:,0], nuz[:,1], nuz[:,2], nuz[:,3]
    px_i, px_ii, px_iii, px_iv = px[:,0], px[:,1], px[:,2], px[:,3]
    pz_i, pz_ii, pz_iii, pz_iv = pz[:,0], pz[:,1], pz[:,2], pz[:,3]
    px_0_i, px_0_ii, px_0_iii, px_0_iv = px_0[:,0], px_0[:,1], px_0[:,2], px_0[:,3]
    pz_0_i, pz_0_ii, pz_0_iii, pz_0_iv = pz_0[:,0], pz_0[:,1], pz_0[:,2], pz_0[:,3]


T1_i, T1_iv =  nT_1[-1,0] - nT_1_del[-1,0], nT_1[-1,3] - nT_1_del[-1,3]
T2_i, T2_iv =  nT_1_del[-1,0], nT_1_del[-1,3]
T3_i, T3_iv =  nT_2[-1,0] - nT_2_del[-1,0], nT_2[-1,3] - nT_2_del[-1,3]
T4_i, T4_iv =  nT_2_del[-1,0], nT_2_del[-1,3]

print (pz_i[-1]-pz_0_i[-1],pz_iv[-1]-pz_0_iv[-1])





if pl.type == 'T_ts':
    plt.figure(figsize=(5,3))
    plt.plot(op.t, nT_1_iv - nT_1_del_iv, label = r'$(\vec{u} \cdot \nabla)T_0$', color = 'k')
    plt.plot(op.t, nT_1_del_iv, label = r'$(\vec{u} \cdot \nabla) \delta T$', color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t, -nT_2_iv + nT_2_del_iv, label = r'$-(\gamma-1)(\nabla \cdot \vec{u})T_0$', color = 'r')
    plt.plot(op.t, -nT_2_del_iv, label = r'$-(\gamma-1)(\nabla \cdot \vec{u}) \delta T$', color = 'r', linestyle = 'dashed', dashes=(3, 2))
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    plt.ylim(-4e-2,1e-2)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc='center right', ncol=2, prop={'size': 12}, frameon=False) 
    plt.title(r'$\mathrm{At} \ (3/4,3/4)$', pad = 10)
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/T_ts_' + pl.point + '.png', dpi =200)
    pass 

elif pl.type == 'ux_ts':
    plt.figure(figsize=(5,3))
    plt.plot(op.t[::100], nux_iv[::100], label = r'$(\vec{u} \cdot \nabla)u_x$', color = 'k')
    plt.plot(op.t[::100], px_0_iv[::100], label = r'$-\partial_x p_0 / \rho$', color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], px_iv[::100] - px_0_iv[::100], label = r'$-\partial_x \delta p / \rho$', color = 'r')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    # plt.ylim(-1.5e-2,1.5e-2)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc='center right', ncol=2, prop={'size': 12}, frameon=False) 
    plt.title(r'$\mathrm{At} \ (3/4,3/4)$', pad = 10)
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/ux_ts_' + pl.point + '.png', dpi =200)
    pass 

elif pl.type == 'uz_ts':
    plt.figure(figsize=(5,3))
    plt.plot(op.t[::100], nuz_iii[::100], label = r'$(\vec{u} \cdot \nabla)u_z$', color = 'k')
    plt.plot(op.t[::100], pz_0_iii[::100] - grid.C3, label = r'$-\partial_z p_0 / \rho - g$', color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], pz_iii[::100] - pz_0_iii[::100], label = r'$-\partial_z \delta p / \rho $', color = 'r')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    # plt.ylim(-8.5e-3,8e-3)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(loc='center right', ncol=2, prop={'size': 12}, frameon=False) 
    plt.title(r'$\mathrm{At} \ (3/4,1/4)$', pad = 10)
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/uz_ts_' + pl.point + '.png', dpi =200)
    pass 

elif pl.type == 'T_ts_old':
    plt.figure(figsize=(5,3))
    line1, = plt.plot(op.t, nT_1_i, label = r'$(1/4,1/4)$', color = 'k')
    line2, = plt.plot(op.t, nT_1_ii, label = r'$(1/4,3/4)$', color = 'r')
    line3, = plt.plot(op.t, nT_1_iii, label = r'$(3/4,1/4)$',color = 'b')
    line4, = plt.plot(op.t, nT_1_iv, label = r'$(3/4,3/4)$',color = 'g')
    plt.plot(op.t, -nT_2_i, color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t, -nT_2_ii, color = 'r', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t, -nT_2_iii, color = 'b', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t, -nT_2_iv, color = 'g', linestyle = 'dashed', dashes=(3, 2))
    line5, = plt.plot([],[], color = 'k', label = r'$(\vec{u} \cdot \nabla)T$')
    line6, = plt.plot([],[], color = 'k', linestyle = '--',dashes=(3, 2), label = r'$-(\gamma-1)(\nabla \cdot \vec{u})T$')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    # plt.ylim(0,1)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    first_legend = plt.legend(handles=[line1,line2,line3,line4],loc='center right', prop={'size': 11}, frameon=False) 
    plt.gca().add_artist(first_legend)
    plt.legend(handles=[line5,line6],loc='center', prop={'size': 12}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/T_ts_old.png', dpi =200)
    pass 

elif pl.type == 'ux_ts_old':
    plt.figure(figsize=(5,3))
    line1, = plt.plot(op.t[::100], nux_i[::100], label = r'$(1/4,1/4)$', color = 'k')
    line2, = plt.plot(op.t[::100], nux_ii[::100], label = r'$(1/4,3/4)$', color = 'r')
    line3, = plt.plot(op.t[::100], nux_iii[::100], label = r'$(3/4,1/4)$',color = 'b')
    line4, = plt.plot(op.t[::100], nux_iv[::100], label = r'$(3/4,3/4)$',color = 'g')
    plt.plot(op.t[::100], px_i[::100], color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], px_ii[::100], color = 'r', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], px_iii[::100], color = 'b', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], px_iv[::100], color = 'g', linestyle = 'dashed', dashes=(3, 2))
    line5, = plt.plot([],[], color = 'k', label = r'$(\vec{u} \cdot \nabla)u_x$')
    line6, = plt.plot([],[], color = 'k', linestyle = '--',dashes=(3, 2), label = r'$-\partial_x p / \rho$')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    plt.ylim(-6e-3,8e-3)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    first_legend = plt.legend(handles=[line1,line2,line3,line4], ncol =1, loc='center right', prop={'size': 11}, frameon=False) 
    plt.gca().add_artist(first_legend)
    plt.legend(handles=[line5,line6],loc='center', prop={'size': 12}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/ux_ts_old.png', dpi =200)
    pass 

elif pl.type == 'uz_ts_old':
    plt.figure(figsize=(5,3))
    line1, = plt.plot(op.t[::100], nuz_i[::100], label = r'$(1/4,1/4)$', color = 'k')
    line2, = plt.plot(op.t[::100], nuz_ii[::100], label = r'$(1/4,3/4)$', color = 'r')
    line3, = plt.plot(op.t[::100], nuz_iii[::100], label = r'$(3/4,1/4)$',color = 'b')
    line4, = plt.plot(op.t[::100], nuz_iv[::100], label = r'$(3/4,3/4)$',color = 'g')
    plt.plot(op.t[::100], pz_i[::100] - 1, color = 'k', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], pz_ii[::100] - 1, color = 'r', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], pz_iii[::100] - 1, color = 'b', linestyle = 'dashed', dashes=(3, 2))
    plt.plot(op.t[::100], pz_iv[::100] - 1, color = 'g', linestyle = 'dashed', dashes=(3, 2))
    line5, = plt.plot([],[], color = 'k', label = r'$(\vec{u} \cdot \nabla)u_z$')
    line6, = plt.plot([],[], color = 'k', linestyle = '--',dashes=(3, 2), label = r'$-\partial_z p / \rho - g$')
    plt.xlabel(r'$t$', fontsize=13)        
    # plt.ylabel(r'$$', fontsize=13) 
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))    
    plt.xlim(0, para.tfinal)      
    # plt.ylim(-6e-3,8e-3)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    first_legend = plt.legend(handles=[line1,line2,line3,line4], ncol =1, loc='center right', prop={'size': 11}, frameon=False) 
    plt.gca().add_artist(first_legend)
    plt.legend(handles=[line5,line6],loc='center', prop={'size': 12}, frameon=False) 
    plt.tight_layout()
    plt.rc('text', usetex=True)          
    plt.savefig(pl.path + pl.output_folder + 'plots/uz_ts_old.png', dpi =200)
    pass 


