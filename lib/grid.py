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
if dim == 2:
    Nx = A*Nz                                   # Number of grid points in x-direction

    Lx = A                                      # Length of box in x-direction

    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
    dx = Lx/Nx                                  # Length between two consecutive grid points in x-direction

    Z = ncp.arange(0,Nz+1)*dz                   # Array consisting of points in z-direction
    X = ncp.arange(0,Nx+A)*dx                   # Array consisting of points in x-direction

    X_mesh, Z_mesh = ncp.meshgrid(X, Z,indexing = 'ij')       # Meshgrids

    temp_d = ncp.zeros_like(Z_mesh)             # Temporary array for derivatives
    pass

elif dim == 1:
    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction         
    Z_mesh = ncp.arange(0,Nz+1)*dz              # Array consisting of points in z-direction
    temp_d = ncp.zeros_like(Z_mesh)             # Temporary array for derivatives
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

alpha = 1/(gamma - 1)                        # Adiabatic index

if reference == 'b':
    # For bottom boundary reference
    T_0 = 1 - epsilon*Z_mesh                 # Temperature constant gradient profile
    rho_0 = T_0**m                           # Density in equilibrium profile
    pass

elif reference == 't':
    # For top boundary reference
    T_0 = 1 + epsilon*(1-Z_mesh)             # Temperature constant gradient profile
    rho_0 = T_0**m                           # Density in equilibrium profile

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
    pass

elif reference == 't':
    t_nu = ncp.sqrt(Ra*(1+alpha)/(Pr*(alpha-m)))  # Viscous time scale
    t_kappa = ncp.sqrt(Pr*Ra*(1+alpha)/(alpha-m)) # Thermal diffusive time scale
    Chi = (1+epsilon)**m                          # Initial density constrast (rho_b/rho_t)
    Ma = epsilon*ncp.sqrt((m+1)/gamma)            # Mach number, U/sqrt(gamma*R_gas*T_t) ; U = sqrt(epsilon*g*l)
    pass

# -------------------------------------------------------------------------- #

# Printing parameters in output
print('# dim =', dim, 'order =', order, 'reference =', reference, 'boundary_u =', boundary_u, 'boundary_T =', boundary_T, 'Scheme =', Scheme)
print('# Nz =', Nz, 'dt =', dt, 'A =', A, 'epsilon =', epsilon, 'alpha =', alpha, 'm =', m, 'Pr =', Pr, 'Ra =', Ra)
print('# t_nu =', ncp.round(t_nu,2), 't_kappa =', ncp.round(t_kappa,2), 'Chi =',  ncp.round(Chi,2), 'Ma =', ncp.round(Ma,2), 'C1 =', C1, 'C2 =', C2, 'C3 =', C3, 'C4 =', C4)