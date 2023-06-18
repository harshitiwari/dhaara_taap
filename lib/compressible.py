from grid import ncp,copy
import grid
from derivative import *
import h5py

class Compressible:

    def __init__(self):
        
        # Velocity vector components
        self.ux = []
        self.uy = []
        self.uz = []

        # Scalar fields
        self.rho = []
        self.T = []

        # Conserved variables
        self.Q = []
        self.F = []
        self.G = []
        self.H = []
        self.E = []

        # For RK time-stepping
        self.Q_copy1 = []
        self.Q_copy2 = []

        # Temporary array
        self.temp = []

        pass

    def set_arrays(self):
        # Setting arrays initially in the memory

        if para.dim == 3:
            # Velocity vector components
            self.ux = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.uy = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.uz = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])

            # Scalar fields
            self.rho = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.T = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])

            # Conserved variables
            self.Q = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.F = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.G = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            self.H = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])

            if para.Scheme == 'RK2':
                self.Q_copy1 = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
                pass
            
            elif para.Scheme == 'RK3':
                self.Q_copy1 = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
                self.Q_copy2 = ncp.zeros([5,grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])

            self.temp = ncp.zeros([grid.Nx+para.A,grid.Ny+para.B,grid.Nz+1])
            pass

        elif para.dim == 2:
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

            if para.coordinates == 'polar':
                self.E = ncp.zeros([4,grid.Nx+para.A,grid.Nz+1])
                pass

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

    def update_conserved(self):
        # Put the primitive variables inside Q's

        if para.dim == 3:
            self.Q[0] = copy(self.rho)
            self.Q[1] = para.A*self.rho*self.ux
            self.Q[2] = para.B*self.rho*self.uy
            self.Q[3] = self.rho*self.uz
            self.Q[4] = self.rho*((1/2)*(para.A**2*self.ux**2 + para.B**2*self.uy**2 + self.uz**2) + grid.C1/(para.gamma-1)*self.T)
            pass

        elif para.dim == 2:
            self.Q[0] = copy(self.rho)
            self.Q[1] = para.A*self.rho*(self.ux*grid.theta_s)
            self.Q[2] = self.rho*self.uz
            self.Q[3] = self.rho*((1/2)*(para.A**2*(self.ux*grid.theta_s)**2 + self.uz**2) + grid.C1/(para.gamma-1)*self.T)
            pass

        elif para.dim == 1:
            self.Q[0] = copy(self.rho)
            self.Q[1] = self.rho*self.uz
            self.Q[2] = self.rho*((1/2)*(self.uz**2) + grid.C1/(para.gamma-1)*self.T)
            pass

        pass

    def compute_convective_flux(self):
        # Convective flux terms

        if para.dim == 3:
            self.F[0] = para.A*self.rho*self.ux
            self.F[1] = self.rho*(para.A**2*self.ux**2 + grid.C1*self.T)
            self.F[2] = para.A*para.B*self.rho*self.ux*self.uy
            self.F[3] = para.A*self.rho*self.ux*self.uz
            self.F[4] = para.A*self.rho*self.ux*((1/2)*(para.A**2*self.ux**2 + para.B**2*self.uy**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)

            self.G[0] = para.B*self.rho*self.uy
            self.G[1] = para.A*para.B*self.rho*self.ux*self.uy
            self.G[2] = self.rho*(para.B**2*self.uy**2 + grid.C1*self.T)
            self.G[3] = para.B*self.rho*self.uy*self.uz
            self.G[4] = para.B*self.rho*self.uy*((1/2)*(para.A**2*self.ux**2 + para.B**2*self.uy**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)

            self.H[0] = self.rho*self.uz
            self.H[1] = para.A*self.rho*self.ux*self.uz
            self.H[2] = para.B*self.rho*self.uy*self.uz
            self.H[2] = self.rho*(self.uz**2 + grid.C1*self.T)
            self.H[3] = self.rho*self.uz*((1/2)*(para.A**2*self.ux**2 + para.B**2*self.uy**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)
            pass

        elif para.dim == 2:
            self.F[0] = para.A*self.rho*(self.ux*grid.theta_s)
            self.F[1] = self.rho*(para.A**2*(self.ux*grid.theta_s)**2 + grid.C1*self.T)
            self.F[2] = para.A*self.rho*(self.ux*grid.theta_s)*self.uz
            self.F[3] = para.A*self.rho*(self.ux*grid.theta_s)*((1/2)*(para.A**2*(self.ux*grid.theta_s)**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)

            self.H[0] = self.rho*self.uz
            self.H[1] = para.A*self.rho*(self.ux*grid.theta_s)*self.uz
            self.H[2] = self.rho*(self.uz**2 + grid.C1*self.T)
            self.H[3] = self.rho*self.uz*((1/2)*(para.A**2*(self.ux*grid.theta_s)**2 + self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)
            pass
        
        elif para.dim == 1:
            self.H[0] = self.rho*self.uz
            self.H[1] = self.rho*(self.uz**2 + grid.C1*self.T)
            self.H[2] = self.rho*self.uz*((1/2)*(self.uz**2) + (para.gamma/(para.gamma - 1))*grid.C1*self.T)
            pass

        pass

    def compute_viscous_flux(self):
        # Viscous flux terms

        if para.dim == 3:
            # Computing tau_xx, tau_yy and tau_zz terms
            self.temp = copy(dfx_c(self.ux))

            self.F[1] -= (4/3)*grid.C2*self.temp                                     
            self.F[4] -= (4/3)*grid.C2*para.A*self.ux*self.temp  
            self.G[2] += (2/3)*grid.C2*self.temp 
            self.G[4] += (2/3)*grid.C2*para.B*self.uy*self.temp                      
            self.H[3] += (2/3)*grid.C2*self.temp                                       
            self.H[4] += (2/3)*grid.C2*self.uz*self.temp

            self.temp = copy(dfy_c(self.uy))

            self.F[1] += (2/3)*grid.C2*self.temp                                     
            self.F[4] += (2/3)*grid.C2*para.A*self.ux*self.temp  
            self.G[2] -= (4/3)*grid.C2*self.temp 
            self.G[4] -= (4/3)*grid.C2*para.B*self.uy*self.temp                      
            self.H[3] += (2/3)*grid.C2*self.temp                                       
            self.H[4] += (2/3)*grid.C2*self.uz*self.temp

            self.temp = copy(dfz_c(self.uz))

            self.F[1] += (2/3)*grid.C2*self.temp                                     
            self.F[4] += (2/3)*grid.C2*para.A*self.ux*self.temp  
            self.G[2] += (2/3)*grid.C2*self.temp 
            self.G[4] += (2/3)*grid.C2*para.B*self.uy*self.temp                      
            self.H[3] -= (4/3)*grid.C2*self.temp                                       
            self.H[4] -= (4/3)*grid.C2*self.uz*self.temp

            # Computing tau_xy terms
            self.temp = copy(dfx_c(self.uy))

            self.F[2] -= grid.C2*(para.B/para.A)*self.temp
            self.F[4] -= grid.C2*(para.B**2/para.A)*self.uy*self.temp 
            self.G[1] -= grid.C2*(para.B/para.A)*self.temp
            self.G[4] -= grid.C2*para.B*self.ux*self.temp

            self.temp = copy(dfy_c(self.ux))

            self.F[2] -= grid.C2*(para.A/para.B)*self.temp
            self.F[4] -= grid.C2*para.A*self.uy*self.temp 
            self.G[1] -= grid.C2*(para.A/para.B)*self.temp
            self.G[4] -= grid.C2*(para.A**2/para.B)*self.ux*self.temp

            # Computing tau_xz terms
            self.temp = copy(dfx_c(self.uz))

            self.F[3] -= grid.C2*self.temp/para.A
            self.F[4] -= grid.C2*self.uz*self.temp/para.A
            self.H[1] -= grid.C2*self.temp/para.A
            self.H[4] -= grid.C2*self.ux*self.temp

            self.temp = copy(dfz_c(self.ux))

            self.F[3] -= grid.C2*para.A*self.temp
            self.F[4] -= grid.C2*para.A*self.uz*self.temp
            self.H[1] -= grid.C2*para.A*self.temp
            self.H[4] -= grid.C2*para.A**2*self.ux*self.temp

            # Computing tau_yz terms
            self.temp = copy(dfy_c(self.uz))

            self.G[3] -= grid.C2*self.temp/para.B
            self.G[4] -= grid.C2*self.uz*self.temp/para.B
            self.H[2] -= grid.C2*self.temp/para.B
            self.H[4] -= grid.C2*self.uy*self.temp

            self.temp = copy(dfz_c(self.uy))

            self.G[3] -= grid.C2*para.B*self.temp
            self.G[4] -= grid.C2*para.B*self.uz*self.temp
            self.H[2] -= grid.C2*para.B*self.temp
            self.H[4] -= grid.C2*para.B**2*self.uy*self.temp

            # Computing q_x, q_y and q_z
            self.temp = copy(dfx_c(self.T))
            self.F[4] -= grid.C4*self.temp/para.A

            self.temp = copy(dfy_c(self.T))
            self.H[4] -= grid.C4*self.temp/para.B

            self.temp = copy(dfz_c(self.T))
            self.H[4] -= grid.C4*self.temp

            pass
        
        elif para.dim == 2:
            # Computing tau_xx and tau_zz terms
            self.temp = copy(dfx_c((self.ux*grid.theta_s)))

            if para.coordinates == 'cartesian':
                self.F[1] -= (4/3)*grid.C2*self.temp                                     
                self.F[3] -= (4/3)*grid.C2*para.A*(self.ux*grid.theta_s)*self.temp                         
                self.H[2] += (2/3)*grid.C2*self.temp                                       
                self.H[3] += (2/3)*grid.C2*self.uz*self.temp  
                pass
            elif para.coordinates == 'polar':
                self.F[1] -= (4/3)*grid.C2*self.temp/grid.Z_mesh                                     
                self.F[3] -= (4/3)*grid.C2*(self.ux*grid.theta_s)*self.temp/grid.Z_mesh                        
                self.H[2] += (2/3)*grid.C2*self.temp/grid.Z_mesh                                       
                self.H[3] += (2/3)*grid.C2*self.uz*self.temp/grid.Z_mesh  

                self.E[2] = (4/3)*grid.C2*self.temp/grid.Z_mesh
                pass

            self.temp = copy(dfz_c(self.uz))

            self.F[1] += (2/3)*grid.C2*self.temp                                     
            self.F[3] += (2/3)*grid.C2*para.A*(self.ux*grid.theta_s)*self.temp                         
            self.H[2] -= (4/3)*grid.C2*self.temp                                       
            self.H[3] -= (4/3)*grid.C2*self.uz*self.temp

            if para.coordinates == 'polar':
                self.E[2] -= (2/3)*grid.C2*self.temp
                pass

            # Computing tau_xz terms
            self.temp = copy(dfz_c((self.ux*grid.theta_s)))

            self.F[2] -= grid.C2*para.A*self.temp                                     
            self.F[3] -= grid.C2*para.A*self.uz*self.temp                          
            self.H[1] -= grid.C2*para.A*self.temp                                       
            self.H[3] -= grid.C2*para.A**2*(self.ux*grid.theta_s)*self.temp    
            
            self.temp = copy(dfx_c(self.uz))

            if para.coordinates == 'cartesian':
                self.F[2] -= grid.C2*self.temp/para.A                                    
                self.F[3] -= grid.C2*self.uz*self.temp/para.A                          
                self.H[1] -= grid.C2*self.temp/para.A                                     
                self.H[3] -= grid.C2*(self.ux*grid.theta_s)*self.temp
                pass
            elif para.coordinates == 'polar':
                self.F[2] -= grid.C2*self.temp/grid.Z_mesh                                     
                self.F[3] -= grid.C2*self.uz*self.temp/grid.Z_mesh                           
                self.H[1] -= grid.C2*self.temp/grid.Z_mesh                                      
                self.H[3] -= grid.C2*(self.ux*grid.theta_s)*self.temp/grid.Z_mesh 
                pass

            # Computing q_x and q_z
            self.temp = copy(dfx_c(self.T))

            if para.coordinates == 'cartesian':
                self.F[3] -= grid.C4*self.temp/para.A
                pass
            elif para.coordinates == 'polar':
                self.F[3] -= grid.C4*self.temp/grid.Z_mesh
                pass

            self.temp = copy(dfz_c(self.T))
            self.H[3] -= grid.C4*self.temp
            pass

        elif para.dim == 1:
            # Computing tau_zz terms
            self.temp = copy(dfz_c(self.uz))
                        
            self.H[1] -= (4/3)*grid.C2*self.temp                                       
            self.H[2] -= (4/3)*grid.C2*self.uz*self.temp

            # Computing q_z 
            self.temp = copy(dfz_c(self.T))
            self.H[2] -= grid.C4*self.temp
            pass

        pass

    def compute_polar_extra_terms(self):
        # Extra flux terms in polar coordinates

        if para.dim == 2:
            self.F[1] -= (4/3)*grid.C2*self.uz/grid.Z_mesh
            self.F[2] += grid.C2*(self.ux*grid.theta_s)/grid.Z_mesh
            self.F[3] -= grid.C2*((1/3)*(self.ux*grid.theta_s)*self.uz)/grid.Z_mesh

            self.H[1] += grid.C2*(self.ux*grid.theta_s)/grid.Z_mesh
            self.H[2] += (2/3)*grid.C2*self.uz/grid.Z_mesh
            self.H[3] += grid.C2*((2/3)*self.uz**2 + (self.ux*grid.theta_s)**2)/grid.Z_mesh

            self.E[2] += (4/3)*grid.C2*self.uz/grid.Z_mesh
            self.E[0] = self.H[0]/grid.Z_mesh
            self.E[1] = 2*self.H[1]/grid.Z_mesh
            self.E[2] += self.H[2]
            self.E[2] -= self.rho*((self.ux*grid.theta_s)**2 + grid.C1*self.T) 
            self.E[2] /= grid.Z_mesh
            self.E[3] = self.H[3]/grid.Z_mesh
            pass

        pass

    def flux_derivative(self):
        # Derivatives of total flux terms

        if para.dim == 3:
            self.F[0] = copy(dfx_c(self.F[0]))
            self.F[1] = copy(dfx_c(self.F[1]))
            self.F[2] = copy(dfx_c(self.F[2]))
            self.F[3] = copy(dfx_c(self.F[3]))
            self.F[4] = copy(dfx_c(self.F[4]))

            self.G[0] = copy(dfy_c(self.G[0]))
            self.G[1] = copy(dfy_c(self.G[1]))
            self.G[2] = copy(dfy_c(self.G[2]))
            self.G[3] = copy(dfy_c(self.G[3]))
            self.G[4] = copy(dfy_c(self.G[4]))

            self.H[0] = copy(dfz_c(self.H[0]))
            self.H[1] = copy(dfz_c(self.H[1]))
            self.H[2] = copy(dfz_c(self.H[2]))
            self.H[3] = copy(dfz_c(self.H[3]))
            self.H[4] = copy(dfz_c(self.H[4]))

            # Gravity source terms
            self.H[3] += grid.C3*self.rho
            self.H[4] += grid.C3*self.rho*self.uz
            pass

        elif para.dim == 2:
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
        # Computing rhs
        self.compute_convective_flux()
        self.compute_viscous_flux()
        if para.coordinates == 'polar':
            self.compute_polar_extra_terms()
            pass
        self.flux_derivative()
        pass

