from grid import ncp,copy
import grid
from compressible import Compressible
from derivative import *
from boundary import *
from saving import print_output, save_fields, timeseries, timeseries_save
import h5py

import sys
sys.path.append("..")

import para

def update_primitive(compress):
    # Update the primitive variables using Q's

    if para.dim == 3:
        compress.rho = copy(compress.Q[0])

        compress.ux = compress.Q[1]/(compress.rho*para.A)
        compress.uy = compress.Q[2]/(compress.rho*para.B)
        compress.uz = compress.Q[3]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[4]/compress.rho) - ((1/2)*(para.A**2*compress.ux**2 + para.B**2*compress.uy**2 + compress.uz**2)))*(para.gamma-1)/grid.C1
        imposeBC_T(compress)
        pass

    elif para.dim == 2:
        compress.rho = copy(compress.Q[0])

        compress.ux = compress.Q[1]/(compress.rho*para.A*grid.theta_s)
        compress.uz = compress.Q[2]/compress.rho
        imposeBC_u(compress)

        compress.T = ((compress.Q[3]/compress.rho) - ((1/2)*(compress.ux**2*para.A**2*grid.theta_s**2 + compress.uz**2)))*(para.gamma-1)/grid.C1
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

    if para.dim == 3:
        compress.Q += -(compress.F/para.A + compress.G/para.B + compress.H)*dt
        pass

    elif para.dim == 2:
        if para.coordinates == 'cartesian':
            compress.Q += -(compress.F/para.A + compress.H)*dt
            pass
        elif para.coordinates == 'polar':
            compress.Q += -(compress.F/para.A + compress.H + compress.E)*dt
            pass
        pass
    
    elif para.dim == 1:
        compress.Q += -compress.H*dt
        pass

    return compress

def time_advance_euler(compress=Compressible()):
    # Time advance using Euler method, 1st order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            print_output(t,compress)
            timeseries(compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            save_fields(t,compress)
            j=j+1
            pass

        compress.update_conserved()
        compress.compute_rhs()
        compress = time_advance_single_step(para.dt, compress)
        update_primitive(compress)

        t = ncp.round(t+para.dt, grid.n1) 
        pass
    timeseries_save(compress)
    pass

def time_advance_RK2(compress=Compressible()):
    # Time advance using RK2 method, 2nd order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            print_output(t,compress)
            timeseries(compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            save_fields(t,compress)
            j=j+1
            pass

        # Saving Q of first point
        compress.update_conserved()
        compress.Q_copy1 = copy(compress.Q)

        # Finding the Q's of middle point 
        compress.compute_rhs()
        compress=time_advance_single_step(para.dt/2, compress)
        update_primitive(compress)

        # Finding the Q's of next point
        compress.compute_rhs()
        compress.Q = copy(compress.Q_copy1)
        compress=time_advance_single_step(para.dt, compress)
        update_primitive(compress)

        t = ncp.round(t+para.dt, grid.n1)
        pass
    timeseries_save(compress)
    pass

def time_advance_RK3(compress=Compressible()):
    # Time advance using RK3 method, 4th order

    t = para.tinit
    j,k = 0,0

    while t <= para.tfinal:
        if ((grid.t_p_step[k]-t)/para.dt) <= para.dt:
            print_output(t,compress)
            timeseries(compress)
            k=k+1
            pass
        if ((grid.t_f_step[j]-t)/para.dt) <= para.dt:
            save_fields(t,compress)
            j=j+1
            pass

        # Saving Q of first point
        compress.update_conserved()
        compress.Q_copy1 = copy(compress.Q)

        # Finding the Q's of 1/3 point and saving for k1
        compress.compute_rhs()
        compress=time_advance_single_step(para.dt/3, compress)
        update_primitive(compress)
        compress.update_conserved()
        compress.Q_copy2 = copy(compress.Q)

        # Finding the Q's of 2/3 point
        compress.Q = copy(compress.Q_copy1)
        compress.compute_rhs()
        compress=time_advance_single_step((2/3)*para.dt, compress)
        update_primitive(compress)

        # Finding the Q's of next point
        compress.Q = (1/4)*compress.Q_copy1 + (3/4)*compress.Q_copy2
        compress.compute_rhs()
        compress=time_advance_single_step((3/4)*para.dt, compress)
        update_primitive(compress)

        t = ncp.round(t+para.dt, grid.n1)
        pass
    timeseries_save(compress)
    pass

