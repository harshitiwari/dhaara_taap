device = "CPU"            # Use Target to set device as CPU or GPU                                        
device_rank = 3           # Set GPU device rank     
 
output_dir = 'output'     # Output directory

dim = 2                   # Dimension of system, 1 or 2            
order = 2                 # Central difference scheme order 2 or 4 
reference = 't'           # Reference values of temperature and density, t: top, b: bottom

boundary_u = 'PFS'        # Boundary condtions for velocity, FS = Free-slip, NS = No-slip, PFS = Periodic in x and free-slip in z, PNS = Periodic in x and no-slip in z
boundary_T = 'PD'         # Boundary condtions for temperature, AD = Adiabatic, PD = Periodic in x

Scheme = 'RK3'            # Time integration scheme, EULER, RK2, RK3

t_f = 20                  # Field files saving time
t_p = 1e-2                # Energies and other parameters saving time

# ----------------------------- Grid parameters ---------------------------- #

tinit = 0                 # Initial time
tfinal = 500              # Final time
dt = 1e-3                 # Single time step

Nz = 128                  # Number of grid points in z-direction
Lz = 1                    # Length of box in z-direction     

# -------------------------------------------------------------------------- #

# --------------------------- Control parameters --------------------------- #

A = 1                     # Aspect ratio
epsilon = 1               # Normalised layer thickness parameter, Delta T/T_r, T_r is reference value
gamma = 5/3               # gamma = C_p/C_v 
m = 1.4                   # Polytropic index
Pr = 1                    # Prandtl number at reference boundary
Ra = 1e3                  # Rayleigh number at reference boundary

# -------------------------------------------------------------------------- #



