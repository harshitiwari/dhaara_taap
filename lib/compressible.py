from grid import ncp,copy
import grid
from derivative import *
import h5py

class Compressible:

    def __init__(self):
        
        # Velocity vector components
        self.ux = []
        self.uz = []

        # Scalar fields
        self.rho = []
        self.T = []

        # Conserved variables
        self.Q = []
        self.F = []
        self.H = []

        # For RK time-stepping
        self.Q_copy1 = []
        self.Q_copy2 = []

        # Temporary array
        self.temp = []

        pass

    def set_arrays(self):

        if para.dim == 2:
            # Velocity vector components
            self.ux = ncp.zeros([grid.Nx+para.A,grid.Nz+1])
            self.uz = ncp.zeros([grid.Nx+para.A,grid.Nz+1])

            # Scalar fields
            self.rho = ncp.zeros([grid.Nx+para.A,grid.Nz+1])
            self.T = ncp.zeros([grid.Nx+para.A,grid.Nz+1])

            # Conserved variables
            self.Q = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])
            self.F = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])
            self.H = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])

            if para.Scheme == 'RK2':
                self.Q_copy1 = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])
                pass
            
            elif para.Scheme == 'RK3':
                self.Q_copy1 = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])
                self.Q_copy2 = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])

            self.temp = ncp.zeros([grid.Nx+para.A,grid.Nz+1])
            pass

        elif para.dim == 1:
            # Velocity
            self.uz = ncp.zeros([grid.Nz+1])

            # Scalar fields
            self.rho = ncp.zeros([grid.Nz+1])
            self.T = ncp.zeros([grid.Nz+1])

            # Conserved variables
            self.Q = ncp.zeros([3,grid.Nz+1])
            self.H = ncp.zeros([3,grid.Nz+1])

            if para.Scheme == 'RK2':
                self.Q_copy1 = ncp.zeros([3,grid.Nz+1])
                pass

            elif para.Scheme == 'RK3':
                self.Q_copy1 = ncp.zeros([3,grid.Nz+1])
                self.Q_copy2 = ncp.zeros([3,grid.Nz+1])

            self.temp = ncp.zeros([grid.Nz+1])
            pass

        pass

    def init_fields(self):

        if para.dim == 2:
            # with h5py.File("output/2D_%.2f.h5" %(grid.tinit), "r") as f:
            #     # List all groups
            #     # print("Keys: %s" % f.keys())
                
            #     data_T = list(f.keys())[0]
            #     data_rho = list(f.keys())[1]
            #     data_ux = list(f.keys())[2]
            #     data_uz = list(f.keys())[3]

            #     # Get the data
            #     self.T = ncp.array(f[data_T])
            #     self.rho = ncp.array(f[data_rho])
            #     self.ux = ncp.array(f[data_ux])
            #     self.uz = ncp.array(f[data_uz])
            
            self.T = copy(grid.T_0)
            self.rho = copy(grid.rho_0)
            self.ux = 0.001*ncp.sin(grid.X_mesh*ncp.pi)*ncp.cos(grid.Z_mesh*ncp.pi*para.A)
            self.uz = - 0.001*ncp.sin(grid.X_mesh*ncp.pi)*ncp.sin(grid.Z_mesh*ncp.pi)
            pass

        elif para.dim == 1:
            # with h5py.File("output/1D_%.2f.h5" %(grid.tinit), "r") as f:
            #     # List all groups
            #     # print("Keys: %s" % f.keys())
                
            #     data_T = list(f.keys())[0]
            #     data_rho = list(f.keys())[1]
            #     data_uz = list(f.keys())[2]

            #     # Get the data
            #     self.T = ncp.array(f[data_T])
            #     self.rho = ncp.array(f[data_rho])
            #     self.uz = ncp.array(f[data_uz])
            
            self.T = copy(grid.T_0)
            self.rho = copy(grid.rho_0)
            self.uz = - 0.001*ncp.sin(grid.Z_mesh*ncp.pi)
            pass

        pass

    def update_conserved(self):

        if para.dim == 2:
            self.Q[0] = copy(self.rho)
            self.Q[1] = para.A*self.rho*self.ux
            self.Q[2] = self.rho*self.uz
            self.Q[3] = self.rho*((1/2)*(para.A**2*self.ux**2 + self.uz**2) + grid.C1/(para.gamma-1)*self.T)
            pass

        elif para.dim == 1:
            self.Q[0] = copy(self.rho)
            self.Q[1] = self.rho*self.uz
            self.Q[2] = self.rho*((1/2)*(self.uz**2) + grid.C1/(para.gamma-1)*self.T)
            pass

        pass

    def compute_convective_flux(self):
        # Convective flux terms

        if para.dim == 2:
            self.F[0] = para.A*self.rho*self.ux
            self.F[1] = self.rho*(para.A**2*self.ux**2 + grid.C1*self.T)
            self.F[2] = para.A*self.rho*self.ux*self.uz
            self.F[3] = para.A*self.rho*self.ux*((1/2)*(para.A**2*self.ux**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)

            self.H[0] = self.rho*self.uz
            self.H[1] = para.A*self.rho*self.ux*self.uz
            self.H[2] = self.rho*(self.uz**2 + grid.C1*self.T)
            self.H[3] = self.rho*self.uz*((1/2)*(para.A**2*self.ux**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)
            pass
        
        elif para.dim == 1:
            self.H[0] = self.rho*self.uz
            self.H[1] = self.rho*(self.uz**2 + grid.C1*self.T)
            self.H[2] = self.rho*self.uz*((1/2)*(self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)
            pass

        pass

    def compute_viscous_flux(self):
        # Viscous flux terms
        
        if para.dim == 2:
            self.temp = copy(dfx_c(self.ux))

            self.F[1] -= (4/3)*grid.C2*self.temp                                     
            self.F[3] -= (4/3)*grid.C2*para.A*self.ux*self.temp                         
            self.H[2] += (2/3)*grid.C2*self.temp                                       
            self.H[3] += (2/3)*grid.C2*self.uz*self.temp                       

            self.temp = copy(dfz_c(self.uz))

            self.F[1] += (2/3)*grid.C2*self.temp                                     
            self.F[3] += (2/3)*grid.C2*para.A*self.ux*self.temp                         
            self.H[2] -= (4/3)*grid.C2*self.temp                                       
            self.H[3] -= (4/3)*grid.C2*self.uz*self.temp

            self.temp = copy(dfz_c(self.ux))

            self.F[2] -= grid.C2*para.A*self.temp                                     
            self.F[3] -= grid.C2*para.A*self.uz*self.temp                          
            self.H[1] -= grid.C2*para.A*self.temp                                       
            self.H[3] -= grid.C2*para.A**2*self.ux*self.temp    

            self.temp = copy(dfx_c(self.uz))

            self.F[2] -= grid.C2*self.temp/para.A                                    
            self.F[3] -= grid.C2*self.uz*self.temp/para.A                          
            self.H[1] -= grid.C2*self.temp/para.A                                     
            self.H[3] -= grid.C2*self.ux*self.temp

            self.temp = copy(dfx_c(self.T))
            self.F[3] -= grid.C4*self.temp/para.A

            self.temp = copy(dfz_c(self.T))
            self.H[3] -= grid.C4*self.temp
            pass

        elif para.dim == 1:
            self.temp = copy(dfz_c(self.uz))
                        
            self.H[1] -= (4/3)*grid.C2*self.temp                                       
            self.H[2] -= (4/3)*grid.C2*self.uz*self.temp

            self.temp = copy(dfz_c(self.T))
            self.H[2] -= grid.C4*self.temp
            pass

        pass

    def flux_derivative(self):

        if para.dim == 2:
            self.F[0] = copy(dfx_c(self.F[0]))
            self.F[1] = copy(dfx_c(self.F[1]))
            self.F[2] = copy(dfx_c(self.F[2]))
            self.F[3] = copy(dfx_c(self.F[3]))

            self.H[0] = copy(dfz_c(self.H[0]))
            self.H[1] = copy(dfz_c(self.H[1]))
            self.H[2] = copy(dfz_c(self.H[2]))
            self.H[3] = copy(dfz_c(self.H[3]))

            # Gravity source terms
            self.H[2] += grid.C3*self.rho
            self.H[3] += grid.C3*self.rho*self.uz
            pass

        elif para.dim == 1:

            self.H[0] = copy(dfz_c(self.H[0]))
            self.H[1] = copy(dfz_c(self.H[1]))
            self.H[2] = copy(dfz_c(self.H[2]))

            # Gravity source terms
            self.H[1] += grid.C3*self.rho
            self.H[2] += grid.C3*self.rho*self.uz
            pass

        pass

    def compute_rhs(self):
        self.update_conserved()
        self.compute_convective_flux()
        self.compute_viscous_flux()
        self.flux_derivative()
        pass

