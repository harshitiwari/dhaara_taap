device = "CPU"            # Use Target to set device as CPU or GPU                                        
device_rank = 3           # Set GPU device rank     
 
output_dir = 'output/2D'     # Output directory

dim = 2                   # Dimension of system, 1, 2 or 3            
order = 4                 # Central difference scheme order 2 or 4 
reference = 't'           # Reference values of temperature and density, t: top, b: bottom

coordinates = 'cartesian'     # Coordinates of the simulation box, cartesian and polar; for polar, z = r and theta = x

# For polar coordinates, theta direction is always periodic, so we have for u: PFS = Free-slip, PNS = No-slip at top and bottom shell; and for T: PDF = Periodic in theta 

boundary_u = 'FS'         # Boundary condtions for velocity, FS = Free-slip, NS = No-slip, PFS = Periodic horizontally and free-slip in z, PNS = Periodic horizontally and no-slip in z, PD = Periodic in all walls
boundary_T = 'ADF'        # Boundary condtions for temperature, ADF = Adiabatic horizontally and fixed in z, PDF = Periodic horizontally and fixed in z, PD = Periodic in all walls

profile = 'static'        # Initial profiles; static: equilibrium profiles, continue: continue some run at time placed in tinit

Scheme = 'RK3'            # Time integration scheme, EULER, RK2, RK3

t_f = 10                  # Field files saving time
t_p = 1e-2                # Energies and other parameters saving time

# ----------------------------- Grid parameters ---------------------------- #

tinit = 0                 # Initial time
tfinal = 500              # Final time
dt = 1e-3                 # Single time step

Nz = 128                  # Number of grid points in z-direction
Lz = 1                    # Length of box in z-direction   

# For polar coordinates

n_theta = 4               # Nz times the number of grids for theta, keep it close to 2 pi
r_b = 0.7                 # Radius of bottom shell for polar coordinate, less than unit radius, r_t = 1

# -------------------------------------------------------------------------- #

# --------------------------- Control parameters --------------------------- #

A = 1                     # Aspect ratio in x-direction
B = 1                     # Aspect ratio in y-direction

epsilon = 1               # Normalised layer thickness parameter, Delta T/T_r, T_r is reference value
gamma = 5/3               # gamma = C_p/C_v 
m = 1.4                   # Polytropic index
Pr = 1                    # Prandtl number at reference boundary
Ra = 1e3                  # Rayleigh number at reference boundary

# -------------------------------------------------------------------------- #



