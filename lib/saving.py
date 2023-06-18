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
    if para.dim == 3:
        return ncp.sum((compress.rho)*(1/2)*(para.A**2*compress.ux**2 + para.B**2*compress.uy**2 + compress.uz**2))/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return ncp.sum((compress.rho)*(1/2)*(para.A**2*(compress.ux*grid.theta_s)**2 + compress.uz**2))/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return ncp.sum((compress.rho)*(1/2)*(compress.uz**2))/((grid.Nz+1))

def total_internal_energy(compress=Compressible()):
    # Total kinetic energy integral
    if para.dim == 3:
        return (grid.C1/(para.gamma-1))*ncp.sum((compress.rho)*(compress.T))/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return (grid.C1/(para.gamma-1))*ncp.sum((compress.rho)*(compress.T))/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return (grid.C1/(para.gamma-1))*ncp.sum((compress.rho)*(compress.T))/((grid.Nz+1))

def vrms(compress=Compressible()):
    # V_rms
    if para.dim == 3:
        return ncp.sqrt(ncp.sum((para.A**2*compress.ux**2 + para.B**2*compress.uy**2 + compress.uz**2))/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1)))
    
    elif para.dim == 2:
        return ncp.sqrt(ncp.sum((para.A**2*(compress.ux*grid.theta_s)**2 + compress.uz**2))/((grid.Nx+1)*(grid.Nz+1)))
    
    elif para.dim == 1:
        return ncp.sqrt(ncp.sum((para.A**2*compress.ux**2 + compress.uz**2))/((grid.Nz+1)))

def flux_convective(compress=Compressible()):
    # Volume average of convective flux
    if para.dim == 3:
        return - ncp.sum(compress.rho*compress.uz*(para.gamma/(para.gamma - 1))*grid.C1*compress.T)/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return - ncp.sum(compress.rho*compress.uz*(para.gamma/(para.gamma - 1))*grid.C1*compress.T)/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return - ncp.sum(compress.rho*compress.uz*(para.gamma/(para.gamma - 1))*grid.C1*compress.T)/((grid.Nz+1))

def flux_radiative(compress=Compressible()):
    # Volume average of radiative flux
    if para.dim == 3:
        return ncp.sum(grid.C4*dfz_c(compress.T))/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return ncp.sum(grid.C4*dfz_c(compress.T))/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return ncp.sum(grid.C4*dfz_c(compress.T))/((grid.Nz+1))

def flux_kinetic(compress=Compressible()):
    # Volume average of kinetic flux
    if para.dim == 3:
        return - ncp.sum((1/2)*compress.rho*compress.uz*(para.A**2*compress.ux**2 + para.B**2*compress.uy**2 + compress.uz**2))/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return - ncp.sum((1/2)*compress.rho*compress.uz*(para.A**2*(compress.ux*grid.theta_s)**2 +  compress.uz**2))/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return - ncp.sum((1/2)*compress.rho*compress.uz**3)/((grid.Nz+1))
    
def criteria(compress=Compressible()):
    # Schwarzschild criteria difference
    if para.dim == 3:
        return ncp.sum(ncp.absolute(dfz_c(compress.T)) - grid.gradT_ad)/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))
    
    elif para.dim == 2:
        return ncp.sum(ncp.absolute(dfz_c(compress.T)) - grid.gradT_ad)/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return ncp.sum(ncp.absolute(dfz_c(compress.T)) - grid.gradT_ad)/((grid.Nz+1))
    
def pressure_term(compress=Compressible()):
    # Volume average of pressure gradient term in uz equation
    if para.dim == 3:
        return grid.C1*ncp.sum(dfz_c(compress.rho*compress.T)/compress.rho)/((grid.Nx+1)*(grid.Ny+1)*(grid.Nz+1))

    elif para.dim == 2:
        return grid.C1*ncp.sum(dfz_c(compress.rho*compress.T)/compress.rho)/((grid.Nx+1)*(grid.Nz+1))
    
    elif para.dim == 1:
        return grid.C1*ncp.sum(dfz_c(compress.rho*compress.T)/compress.rho)/((grid.Nz+1))

def nlinux(nx,nz,compress=Compressible()):
    # Non-linear term in ux equation
    if para.dim == 2:
        return (grid.theta_s*compress.ux[nx,nz]*dfx_p(grid.theta_s*compress.ux,nx,nz) + compress.uz[nx,nz]*dfz_p(grid.theta_s*compress.ux,nx,nz))

def nlinuz(nx,nz,compress=Compressible()):
    # Non-linear term in uz equation
    if para.dim == 2:
        return (grid.theta_s*compress.ux[nx,nz]*dfx_p(compress.uz,nx,nz) + compress.uz[nx,nz]*dfz_p(compress.uz,nx,nz))

def gradpx(nx,nz,compress=Compressible()):
    # Pressure gradient term in ux equation
    if para.dim == 2:
        return -(grid.C1/compress.rho[nx,nz])*dfx_p(compress.rho*compress.T,nx,nz)

def gradpz(nx,nz,compress=Compressible()):
    # Pressure gradient term in uz equation
    if para.dim == 2:
        return -(grid.C1/compress.rho[nx,nz])*dfz_p(compress.rho*compress.T,nx,nz)

def gradpx_0(nx,nz,compress=Compressible()):
    # Pressure gradient term in ux equation
    if para.dim == 2:
        return -(grid.C1/compress.rho[nx,nz])*dfx_p(grid.rho_0*grid.T_0,nx,nz)

def gradpz_0(nx,nz,compress=Compressible()):
    # Pressure gradient term in uz equation
    if para.dim == 2:
        return -(grid.C1/compress.rho[nx,nz])*dfz_p(grid.rho_0*grid.T_0,nx,nz)

def nlinT_1(nx,nz,compress=Compressible()):
    # Non-linear 1st term in temperature equation at point (nx,nz)
    if para.dim == 2:
        return (grid.theta_s*compress.ux[nx,nz]*dfx_p(compress.T,nx,nz) + compress.uz[nx,nz]*dfz_p(compress.T,nx,nz)) 

def nlinT_2(nx,nz,compress=Compressible()):
    # Non-linear 2nd term in temperature equation at point (nx,nz)
    if para.dim == 2:
        return (para.gamma-1)*compress.T[nx,nz]*(dfx_p(grid.theta_s*compress.ux,nx,nz) + dfz_p(compress.uz,nx,nz)) 

def nlinT_1_del(nx,nz,compress=Compressible()):
    # Non-linear 1st term of temperature perturbation at point (nx,nz)
    if para.dim == 2:
        return (grid.theta_s*compress.ux[nx,nz]*dfx_p(compress.T - grid.T_0,nx,nz) + compress.uz[nx,nz]*dfz_p(compress.T - grid.T_0,nx,nz)) 

def nlinT_2_del(nx,nz,compress=Compressible()):
    # Non-linear 2nd term of temperature perturbation at point (nx,nz)
    if para.dim == 2:
        return (para.gamma-1)*(compress.T[nx,nz] - grid.T_0[nx,nz])*(dfx_p(grid.theta_s*compress.ux,nx,nz) + dfz_p(compress.uz,nx,nz)) 

def print_output(t, compress=Compressible()):

    Ke = total_kinetic_energy(compress)
    Ie = total_internal_energy(compress)
    V_rms = vrms(compress)
    f_c = flux_convective(compress) 
    f_r = flux_radiative(compress) 
    f_k = flux_kinetic(compress)
    cri = criteria(compress)
    p_t = pressure_term(compress)

    print(t, Ke, Ie, V_rms, f_c, f_r, f_k, cri, p_t)

    if math.isnan(Ke):
        print ('# Try different parameters..Energy became infinity, code blew up at t = ', t)    
        sys.exit(1)

    pass

def timeseries(compress=Compressible()):
    # Time series of various terms of our equations at particular 4 points in the box

    nx1, nz1 = int((1/4)*grid.Nx), int((1/4)*para.Nz)
    nx2, nz2 = int((1/4)*grid.Nx), int((3/4)*para.Nz)
    nx3, nz3 = int((3/4)*grid.Nx), int((1/4)*para.Nz)
    nx4, nz4 = int((3/4)*grid.Nx), int((3/4)*para.Nz)

    grid.nux.append([nlinux(nx1,nz1,compress), nlinux(nx2,nz2,compress), nlinux(nx3,nz3,compress), nlinux(nx4,nz4,compress)])
    grid.nuz.append([nlinuz(nx1,nz1,compress), nlinuz(nx2,nz2,compress), nlinuz(nx3,nz3,compress), nlinuz(nx4,nz4,compress)])
    grid.px.append([gradpx(nx1,nz1,compress), gradpx(nx2,nz2,compress), gradpx(nx3,nz3,compress), gradpx(nx4,nz4,compress)])
    grid.pz.append([gradpz(nx1,nz1,compress), gradpz(nx2,nz2,compress), gradpz(nx3,nz3,compress), gradpz(nx4,nz4,compress)])
    grid.px_0.append([gradpx_0(nx1,nz1,compress), gradpx_0(nx2,nz2,compress), gradpx_0(nx3,nz3,compress), gradpx_0(nx4,nz4,compress)])
    grid.pz_0.append([gradpz_0(nx1,nz1,compress), gradpz_0(nx2,nz2,compress), gradpz_0(nx3,nz3,compress), gradpz_0(nx4,nz4,compress)])
    grid.nT_1.append([nlinT_1(nx1,nz1,compress), nlinT_1(nx2,nz2,compress), nlinT_1(nx3,nz3,compress), nlinT_1(nx4,nz4,compress)])
    grid.nT_2.append([nlinT_2(nx1,nz1,compress), nlinT_2(nx2,nz2,compress), nlinT_2(nx3,nz3,compress), nlinT_2(nx4,nz4,compress)])
    grid.nT_1_del.append([nlinT_1_del(nx1,nz1,compress), nlinT_1_del(nx2,nz2,compress), nlinT_1_del(nx3,nz3,compress), nlinT_1_del(nx4,nz4,compress)])
    grid.nT_2_del.append([nlinT_2_del(nx1,nz1,compress), nlinT_2_del(nx2,nz2,compress), nlinT_2_del(nx3,nz3,compress), nlinT_2_del(nx4,nz4,compress)])
    
    pass

def save_fields(t,compress=Compressible()):
    # Saving the fields

    parameters = ncp.array([grid.alpha, para.m, para.epsilon, para.Pr, para.Ra, para.dt])
    if para.dim == 3:
        hf1 = h5py.File(para.output_dir + "/3D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0))
            hf1.create_dataset('ux', data=compress.ux)
            hf1.create_dataset('uy', data=compress.uy)
            hf1.create_dataset('uz', data=compress.uz)
            hf1.create_dataset('T', data=(compress.T - grid.T_0))
            hf1.create_dataset('parameters', data=parameters)
            hf1.close()
            pass
        else:
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0).get())
            hf1.create_dataset('ux', data=compress.ux.get())
            hf1.create_dataset('uy', data=compress.uy.get())
            hf1.create_dataset('uz', data=compress.uz.get())
            hf1.create_dataset('T', data=(compress.T - grid.T_0).get())
            hf1.create_dataset('parameters', data=parameters.get())
            hf1.close()
            pass
        pass
    elif para.dim == 2:
        hf1 = h5py.File(para.output_dir + "/2D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0))
            hf1.create_dataset('ux', data=compress.ux)
            hf1.create_dataset('uz', data=compress.uz)
            hf1.create_dataset('T', data=(compress.T - grid.T_0))
            hf1.create_dataset('parameters', data=parameters)
            hf1.close()
            pass
        else:
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0).get())
            hf1.create_dataset('ux', data=compress.ux.get())
            hf1.create_dataset('uz', data=compress.uz.get())
            hf1.create_dataset('T', data=(compress.T - grid.T_0).get())
            hf1.create_dataset('parameters', data=parameters.get())
            hf1.close()
            pass
        pass
    elif para.dim == 1:
        hf1 = h5py.File(para.output_dir + "/1D_%.2f.h5" % (t), 'w')
        if para.device == 'CPU':
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0))
            hf1.create_dataset('uz', data=compress.uz)
            hf1.create_dataset('T', data=(compress.T - grid.T_0))
            hf1.create_dataset('parameters', data=parameters)
            hf1.close()
            pass
        else:
            hf1.create_dataset('rho', data=(compress.rho - grid.rho_0).get())
            hf1.create_dataset('uz', data=compress.uz.get())
            hf1.create_dataset('T', data=(compress.T - grid.T_0).get())
            hf1.create_dataset('parameters', data=parameters).get()
            hf1.close()
            pass
        pass
    pass

def timeseries_save(compress=Compressible()):
    # Saving the timeseries

    hf2 = h5py.File(para.output_dir + "/timeseries.h5", 'w')
    if para.device == 'CPU':
        hf2.create_dataset('nux', data=grid.nux)
        hf2.create_dataset('nuz', data=grid.nuz)
        hf2.create_dataset('px', data=grid.px)
        hf2.create_dataset('pz', data=grid.pz)
        hf2.create_dataset('px_0', data=grid.px_0)
        hf2.create_dataset('pz_0', data=grid.pz_0)
        hf2.create_dataset('nT_1', data=grid.nT_1)
        hf2.create_dataset('nT_2', data=grid.nT_2)
        hf2.create_dataset('nT_1_del', data=grid.nT_1_del)
        hf2.create_dataset('nT_2_del', data=grid.nT_2_del)
        hf2.close()
        pass
    else:
        hf2.create_dataset('nux', data=grid.nux.get())
        hf2.create_dataset('nuz', data=grid.nuz.get())
        hf2.create_dataset('px', data=grid.px.get())
        hf2.create_dataset('pz', data=grid.pz.get())
        hf2.create_dataset('px_0', data=grid.px_0.get())
        hf2.create_dataset('pz_0', data=grid.pz_0.get())
        hf2.create_dataset('nT_1', data=grid.nT_1.get())
        hf2.create_dataset('nT_2', data=grid.nT_2.get())
        hf2.create_dataset('nT_1_del', data=grid.nT_1_del.get())
        hf2.create_dataset('nT_2_del', data=grid.nT_2_del.get())
        hf2.close()
        pass

    pass
