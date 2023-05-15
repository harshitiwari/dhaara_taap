from grid import ncp,copy
import grid
from compressible import Compressible
from derivative import *
from boundary import *
from saving import print_output, save_fields
import h5py

import sys
sys.path.append("..")

import para

def update_primitive(compress):
    # Update the primitive variables using Q's

    if para.dim == 2:
        compress.rho = copy(compress.Q[0])

        compress.ux = compress.Q[1]/(compress.rho*para.A)
        compress.uz = compress.Q[2]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[3]/compress.rho) - ((1/2)*(compress.ux**2*para.A**2 + compress.uz**2)))*(para.gamma-1)/grid.C1
        imposeBC_T(compress)
        pass

    elif para.dim == 1:
        compress.rho = copy(compress.Q[0])

        compress.uz = compress.Q[1]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[2]/compress.rho) - ((1/2)*(compress.uz**2)))*(para.gamma-1)/grid.C1
        imposeBC_T(compress)
        pass

    pass

def time_advance_single_step(dt, compress=Compressible()):
    # Q for single time step dt

    compress.compute_rhs()

    if para.dim == 2:
        compress.Q -= (compress.F/para.A + compress.H)*dt
        pass
    
    elif para.dim == 1:
        compress.Q -= compress.H*dt
        pass

    update_primitive(compress)

    return compress

def time_advance_euler(compress=Compressible()):
    # Time advance using Euler method, 1st order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            compress = print_output(t,compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            compress = save_fields(t,compress)
            j=j+1
            pass

        compress = time_advance_single_step(para.dt, compress)

        t = ncp.round(t+para.dt, grid.n1) 
        pass

    pass

def time_advance_RK2(compress=Compressible()):
    # Time advance using RK2 method, 2nd order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            print_output(t,compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            save_fields(t,compress)
            j=j+1
            pass

        # Saving Q of first point
        compress.Q_copy1 = copy(compress.Q)
        # Finding the Q's of middle point 
        compress=time_advance_single_step(para.dt/2, compress)
        compress.compute_rhs()

        compress.Q = copy(compress.Q_copy1)  
        # Finding the Q's of next point
        compress=time_advance_single_step(para.dt, compress)

        t = ncp.round(t+para.dt, grid.n1)
        pass

    pass

def time_advance_RK3(compress=Compressible()):
    # Time advance using RK3 method, 4th order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            print_output(t,compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            save_fields(t,compress)
            j=j+1
            pass

        # Saving Q of first point
        compress.Q_copy1 = copy(compress.Q)
        # Finding the Q's of 1/3 point 
        compress=time_advance_single_step(para.dt/2, compress)
        # Finding k2 and saving it to Q_copy2
        compress.compute_rhs()
        compress.Q_copy2 = para.dt*copy(compress.Q)

        # Finding the Q's of 2/3 point
        compress.Q = compress.Q_copy1 + (2/3)*compress.Q_copy2
        # Finding k3 and saving it to Q_copy2
        compress.compute_rhs()
        compress.Q_copy2 = para.dt*copy(compress.Q)

        # Finding k1 and adding it to first point
        compress.Q = copy(compress.Q_copy1)
        compress=time_advance_single_step(para.dt/4, compress)
        compress.Q += (3/4)*compress.Q_copy2 

        t = ncp.round(t+para.dt, grid.n1)
        pass
        
    pass

