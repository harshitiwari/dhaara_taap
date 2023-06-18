import h5py
import sys
sys.path.append("..")
import para
sys.path.append("../lib")
from grid import ncp,copy
import grid
from compressible import Compressible

# Initial field profiles 
def init_fields(compress):

    if para.profile == 'static':
        # Equilibrium profiles for temperature and density
        compress.T = copy(grid.T_0)
        compress.rho = copy(grid.rho_0)

        if para.dim == 3:  
            compress.ux = 0.001*ncp.sin(grid.X_mesh*ncp.pi)*ncp.cos(grid.Y_mesh*ncp.pi)*ncp.cos(grid.Z_mesh*ncp.pi)
            compress.uy = 0.001*ncp.cos(grid.X_mesh*ncp.pi)*ncp.sin(grid.Y_mesh*ncp.pi)*ncp.cos(grid.Z_mesh*ncp.pi)
            compress.uz = - 0.001*ncp.cos(grid.X_mesh*ncp.pi)*ncp.cos(grid.Y_mesh*ncp.pi)*ncp.sin(grid.Z_mesh*ncp.pi)
            pass

        elif para.dim == 2:
            if para.coordinates == 'cartesian':
                compress.ux = 0.001*ncp.sin(grid.X_mesh*ncp.pi)*ncp.cos(grid.Z_mesh*ncp.pi)
                compress.uz = - 0.001*ncp.cos(grid.X_mesh*ncp.pi)*ncp.sin(grid.Z_mesh*ncp.pi)
                pass
            elif para.coordinates == 'polar':
                compress.ux = 0.001*ncp.sin(grid.X_mesh)
                pass
            pass

        elif para.dim == 1:
            compress.uz = - 0.001*ncp.sin(grid.Z_mesh*ncp.pi)
            pass

        pass

    elif para.profile == 'continue':
        # Continue previous run from some time, put it in para.tinit
        with h5py.File("output/2D_%.2f.h5" %(grid.tinit), "r") as f:
                # List all groups
                # print("Keys: %s" % f.keys())
                
                data_T = list(f.keys())[0]
                data_rho = list(f.keys())[2]
                data_ux = list(f.keys())[3]
                data_uz = list(f.keys())[4]

                # Get the data
                compress.T = ncp.array(f[data_T])
                compress.rho = ncp.array(f[data_rho])
                compress.ux = ncp.array(f[data_ux])
                compress.uz = ncp.array(f[data_uz])
        pass

    pass


