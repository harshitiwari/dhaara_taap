from grid import ncp,copy
import grid
from compressible import Compressible
from derivative import *
import h5py
import math as math

import sys
sys.path.append("..")

import para

def total_kinetic_energy(compress=Compressible()):
    # Total kinetic energy integral
    if para.dim == 2:
        return ncp.sum((compress.rho)*(1/2)*(para.A**2*compress.ux**2 + compress.uz**2))/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return ncp.sum((compress.rho)*(1/2)*(compress.uz**2))/(grid.Nz)

def total_internal_energy(compress=Compressible()):
    # Total kinetic energy integral
    if para.dim == 2:
        return (grid.C1/(para.gamma-1))*ncp.sum((compress.rho)*(compress.T))/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return (grid.C1/(para.gamma-1))*ncp.sum((compress.rho)*(compress.T))/(grid.Nz)

def vrms(compress=Compressible()):
    # V_rms
    if para.dim == 2:
        return ncp.sqrt(ncp.sum((para.A**2*compress.ux**2 + compress.uz**2)))/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return ncp.sqrt(ncp.sum((para.A**2*compress.ux**2 + compress.uz**2)))/(grid.Nz)

def Nusselt_number_b(compress=Compressible()):

    if para.dim == 2:
        return - ncp.sum(dfz_c((compress.T))[:,0])/(grid.Nx*para.epsilon)
    
    elif para.dim == 1:
        return - ncp.sum(dfz_c((compress.T))[0])/para.epsilon
    
def Nusselt_number_t(compress=Compressible()):

    if para.dim == 2:
        return - ncp.sum(dfz_c((compress.T))[:,-1])/(grid.Nx*para.epsilon)
    
    elif para.dim == 1:
        return - ncp.sum(dfz_c((compress.T))[-1])/para.epsilon

def criteria(compress=Compressible()):

    if para.dim == 2:
        return ncp.sum(dfz_c(compress.T) + grid.gradT_ad)/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return ncp.sum(dfz_c(compress.T) + grid.gradT_ad)/(grid.Nz)
    
def pressure_term(compress=Compressible()):

    if para.dim == 2:
        return grid.C1*ncp.sum(dfz_c(compress.rho*compress.T))/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return grid.C1*ncp.sum(dfz_c(compress.rho*compress.T))/(grid.Nz)
    
def gravity_term(compress=Compressible()):

    if para.dim == 2:
        return grid.C3*ncp.sum(compress.rho)/(grid.Nx*grid.Nz)
    
    elif para.dim == 1:
        return grid.C3*ncp.sum(compress.rho)/(grid.Nz)

def print_output(t, compress=Compressible()):

    Ke = total_kinetic_energy(compress)
    Ie = total_internal_energy(compress)
    V_rms = vrms(compress)
    Nu_b = Nusselt_number_b(compress) 
    Nu_t = Nusselt_number_t(compress) 
    cri = criteria(compress)
    p_t = pressure_term(compress)
    g_t = gravity_term(compress)

    print(t, Ke, Ie, V_rms, Nu_b, Nu_t, cri, p_t, g_t)

    if math.isnan(Ke):
        print ('# Try different parameters..Energy became infinity, code blew up at t = ', t)    
        sys.exit(1)

    pass

def save_fields(t,compress=Compressible()):

    parameters = ncp.array([grid.alpha, para.m, para.Pr, para.Ra, para.dt])

    if para.dim == 2:
        hf = h5py.File(para.output_dir + "/2D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf.create_dataset('rho', data=(compress.rho - grid.rho_0))
            hf.create_dataset('ux', data=compress.ux)
            hf.create_dataset('uz', data=compress.uz)
            hf.create_dataset('T', data=(compress.T - grid.T_0))
            hf.create_dataset('parameters', data=parameters)
            hf.close()
            pass
        else:
            hf.create_dataset('rho', data=(compress.rho - grid.rho_0).get())
            hf.create_dataset('ux', data=compress.ux.get())
            hf.create_dataset('uz', data=compress.uz.get())
            hf.create_dataset('T', data=(compress.T - grid.T_0).get())
            hf.create_dataset('parameters', data=parameters.get())
            hf.close()
            pass
        pass
    elif para.dim == 1:
        hf = h5py.File(para.output_dir + "/1D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf.create_dataset('rho', data=(compress.rho - grid.rho_0))
            hf.create_dataset('uz', data=compress.uz)
            hf.create_dataset('T', data=(compress.T - grid.T_0))
            hf.create_dataset('parameters', data=parameters)
            hf.close()
            pass
        else:
            hf.create_dataset('rho', data=(compress.rho - grid.rho_0).get())
            hf.create_dataset('uz', data=compress.uz.get())
            hf.create_dataset('T', data=(compress.T - grid.T_0).get())
            hf.create_dataset('parameters', data=parameters).get()
            hf.close()
            pass
        pass
    pass

