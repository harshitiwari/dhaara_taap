import sys
sys.path.append("..")
from para import *

if device == 'CPU':
    import numpy as ncp
    from numpy import copy
    
else:
    import cupy as ncp
    from cupy import copy

    dev = ncp.cuda.Device(device_rank)
    dev.use()

# ----------------------------- Time variables ----------------------------- #

t = ncp.arange(tinit,tfinal+dt,dt)                  # Time axis

Nf = int((tfinal-tinit)/t_f)                        # Number of field saving times
t_f_step = ncp.linspace(tinit,tfinal,Nf+1)          # Field saving times

Np = int((tfinal-tinit)/t_p)                        # Number of field saving times
t_p_step = ncp.linspace(tinit,tfinal,Np+1)          # Field saving times

n1=1                                                # Used for rounding time in time-stepping
n2=10 
while dt*n2!=1:
    n1=n1+1
    n2=n2*10
    
# -------------------------------------------------------------------------- #

# ----------------------------- Grid variables ----------------------------- #

if dim == 3:
    Nx = A*Nz                                   # Number of grid points in x-direction
    Ny = B*Nz                                   # Number of grid points in y-direction
    
    Lx = A                                      # Length of box in x-direction
    Ly = B                                      # Length of box in y-direction

    dx = Lx/Nx                                  # Length between two consecutive grid points in x-direction
    dy = Ly/Ny                                  # Length between two consecutive grid points in y-direction
    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
    
    X = ncp.arange(0,Nx+A)*dx                   # Array consisting of points in x-direction
    Y = ncp.arange(0,Ny+B)*dx                   # Array consisting of points in y-direction
    Z = ncp.arange(0,Nz+1)*dz                   # Array consisting of points in z-direction
    
    X_mesh, Y_mesh, Z_mesh = ncp.meshgrid(X, Y, Z,indexing = 'ij')       # Meshgrids

    temp_dx = ncp.zeros_like(X_mesh)            # Temporary array for x-derivative
    temp_dy = ncp.zeros_like(Y_mesh)            # Temporary array for y-derivative
    temp_dz = ncp.zeros_like(Z_mesh)            # Temporary array for z-derivative

    pass

elif dim == 2:
    if coordinates == 'cartesian':
        # For Cartesian coordinates
        Nx = A*Nz                                   # Number of grid points in x-direction

        Lx = A                                      # Length of box in x-direction

        dx = Lx/Nx                                  # Length between two consecutive grid points in x-direction
        dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
        
        X = ncp.arange(0,Nx+A)*dx                   # Array consisting of points in x-direction
        Z = ncp.arange(0,Nz+1)*dz                   # Array consisting of points in z-direction

        theta_s = 1                                 # Theta length scale, for cartesian its 1
        pass

    elif coordinates == 'polar':
        # For polar coordinates
        Nx = int(n_theta*Nz)                        # Number of grid points in theta-direction

        dz = (1-r_b)/Nz                             # Length between two consecutive grid points in r-direction         
        dx = 2*ncp.pi/Nx                            # Length between two consecutive grid points in theta-direction

        Z = ncp.linspace(r_b,1,Nz+1)                # Array consisting of points in r-direction from r_b to 1
        X = ncp.arange(0,Nx+1)*dx                   # Array consisting of points in theta-direction

        theta_s = 2*ncp.pi                          # Theta length scale, for polar its 2 pi
        pass
    
    X_mesh, Z_mesh = ncp.meshgrid(X, Z,indexing = 'ij')       # Meshgrids

    temp_dx = ncp.zeros_like(X_mesh)                # Temporary array for x-derivative
    temp_dz = ncp.zeros_like(Z_mesh)                # Temporary array for z-derivative

    pass

elif dim == 1:
    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
    Z_mesh = ncp.arange(0,Nz+1)*dz              # Array consisting of points in z-direction
    temp_dz = ncp.zeros_like(Z_mesh)            # Temporary array for z-derivative
    pass

# -------------------------------------------------------------------------- #

# --------------------------- Boundary parameters -------------------------- #

if reference == 'b':
    T_t = 1 - epsilon                        # Temperature at top boundary
    T_b = 1                                  # Temperature at bottom boundary
    pass

elif reference == 't':
    T_t = 1                                  # Temperature at top boundary
    T_b = 1 + epsilon                        # Temperature at bottom boundary
    pass

# -------------------------------------------------------------------------- #

# -------------------------- Equilibrium profile --------------------------- #

alpha = ncp.round(1/(gamma - 1),1)           # Adiabatic index

if reference == 'b':
    # For bottom boundary reference
    T_0 = 1 - epsilon*Z_mesh                 # Temperature constant gradient profile
    rho_0 = T_0**m                           # Density in equilibrium profile
    gradT_ad = epsilon*(1+m)/(1+alpha)       # Adiabatic temperature gradient, g/Cp
    pass

elif reference == 't':
    # For top boundary reference
    if coordinates == 'cartesian':
        T_0 = 1 + epsilon*(1-Z_mesh)         # Temperature constant gradient profile
        rho_0 = T_0**m                       # Density in equilibrium profile
        pass
    elif coordinates == 'polar':
        T_0 = 1 + epsilon*(1-Z_mesh)/(1-r_b) # Temperature constant gradient profile 
        rho_0 = T_0**(m*(1-r_b)-r_b)         # Density in equilibrium profile 
        pass

    gradT_ad = epsilon*(1+m)/(1+alpha)       # Adiabatic temperature gradient, g/Cp
    
    pass

# -------------------------------------------------------------------------- #

# ------------------------- Constants in equations ------------------------- #

if reference == 'b':
    # For lower boundary reference
    C1 = 1/(epsilon**2*(m+1))
    C2 = ncp.sqrt(Pr/Ra)
    C3 = 1/epsilon
    C4 = gamma/((gamma-1)*epsilon**2*(m+1)*ncp.sqrt(Ra*Pr))
    pass

elif reference == 't':
    # For upper boundary reference
    C1 = 1/(epsilon**2*(m+1))
    C2 = ncp.sqrt(Pr*(alpha-m)/(Ra*(1+alpha)))
    C3 = 1/epsilon
    C4 = ncp.sqrt((alpha-m)*(1+alpha)/(Pr*Ra*(1+m)**2*epsilon**4))
    pass

# -------------------------------------------------------------------------- #

# --------------------------- Related parameters --------------------------- #

if reference == 'b':
    t_nu = ncp.sqrt(Ra/Pr)                        # Viscous time scale
    t_kappa = ncp.sqrt(Pr*Ra)                     # Thermal diffusive time scale
    Chi = 1/(1-epsilon)**m                        # Initial density constrast (rho_b/rho_t)
    Ma = epsilon*ncp.sqrt((m+1)/gamma)            # Mach number, U/sqrt(gamma*R_gas*T_b) ; U = sqrt(epsilon*g*l)
    gradT_diff = epsilon*(alpha-m)/(alpha+1)      # Delta T/Lz - g/Cp, Schwarschild criteria temperature gradient difference
    pass

elif reference == 't':
    if coordinates == 'cartesian':
        t_nu = ncp.sqrt(Ra*(1+alpha)/(Pr*(alpha-m)))  # Viscous time scale
        t_kappa = ncp.sqrt(Pr*Ra*(1+alpha)/(alpha-m)) # Thermal diffusive time scale
        Chi = (1+epsilon)**m                          # Initial density constrast (rho_b/rho_t)
        Ma = epsilon*ncp.sqrt((m+1)/gamma)            # Mach number, U/sqrt(gamma*R_gas*T_t) ; U = sqrt(epsilon*g*l)
        gradT_diff = epsilon*(alpha-m)/(alpha+1)      # Delta T/Lz - g/Cp, Schwarschild criteria temperature gradient difference
        pass
    elif coordinates == 'polar':
        t_nu = ncp.sqrt(Ra*(1+alpha)/(Pr*(alpha-m)))                        # Viscous time scale
        t_kappa = ncp.sqrt(Pr*Ra*(1+alpha)/(alpha-m))                       # Thermal diffusive time scale
        Chi = (1+epsilon)**(m*(1-r_b)-r_b)                                  # Initial density constrast (rho_b/rho_t)
        Ma = epsilon*ncp.sqrt((m+1)/gamma)                                  # Mach number, U/sqrt(gamma*R_gas*T_t) ; U = sqrt(epsilon*g*l)
        gradT_diff = epsilon*(alpha-m + r_b*(1-m))/((alpha+1)*(1-r_b))      # Delta T/Lz - g/Cp, Schwarschild criteria temperature gradient difference
        pass
    pass

# -------------------------------------------------------------------------- #

nux, nuz, px, pz, px_0, pz_0, nT_1, nT_2, nT_1_del, nT_2_del = [], [], [], [], [], [], [], [], [], []           # Time series saving list for various terms in our equation

# Printing parameters in output
print('# dim =', dim, ', order =', order, ', reference =', reference, ', coordinates =', coordinates, ', boundary_u =', boundary_u, ', boundary_T =', boundary_T, ', Scheme =', Scheme)
print('\n# Nz =', Nz, ', n_theta =', n_theta, ', r_b =', r_b, ', dt =', dt, ', A =', A, ', B =', B, ', epsilon =', epsilon, ', alpha =', alpha, ', m =', m, ', Pr =', Pr, ', Ra =', Ra)
print('\n# t_nu =', ncp.round(t_nu,2), ', t_kappa =', ncp.round(t_kappa,2), ', Chi =',  ncp.round(Chi,2), ', Ma =', ncp.round(Ma,2),'gradT_diff =', ncp.round(gradT_diff,2), ', C1 =', C1, ', C2 =', C2, ', C3 =', C3, ', C4 =', C4)
print('\n \n# The following columns contains in order: t, Ke, Ie, V_rms, f_c, f_r, f_k, criteria diff, p_av')