import grid
from compressible import Compressible

import sys
sys.path.append("..")

import para

def imposeBC_u(compress = Compressible()):
    # Boundary condition for velocity fields

    if para.dim == 3:
        if para.boundary_u == 'FS':
            # Free-slip boundary condition for velocity at all walls
            compress.ux[0,:,:] = compress.ux[-1,:,:] = 0
            compress.uy[0,:,:] = compress.uy[1,:,:] 
            compress.uy[-1,:,:]  = compress.uy[-2,:,:]
            compress.uz[0,:,:] = compress.uz[1,:,:] 
            compress.uz[-1,:,:]  = compress.uz[-2,:,:]

            compress.ux[:,0,:] = compress.ux[:,1,:] 
            compress.ux[:,-1,:]  = compress.ux[:,-2,:]
            compress.uy[:,0,:] = compress.uy[:,-1,:] = 0
            compress.uz[:,0,:] = compress.uz[:,1,:] 
            compress.uz[:,-1,:]  = compress.uz[:,-2,:]
            
            compress.ux[:,:,0] = compress.ux[:,:,1] 
            compress.ux[:,:,-1] = compress.ux[:,:,-2]
            compress.uy[:,:,0] = compress.uy[:,:,1] 
            compress.uy[:,:,-1] = compress.uy[:,:,-2]       
            compress.uz[:,:,0] = compress.uz[:,:,-1] = 0
            pass
        
        elif para.boundary_u == 'NS':
            # No-slip boundary condition for velocity at all walls
            compress.ux[0,:,:] = compress.ux[-1,:,:] = 0
            compress.uy[0,:,:] = compress.uy[-1,:,:] = 0
            compress.uz[0,:,:] = compress.uz[-1,:,:] = 0

            compress.ux[:,0,:] = compress.ux[:,-1,:] = 0
            compress.uy[:,0,:] = compress.uy[:,-1,:] = 0
            compress.uz[:,0,:] = compress.uz[:,-1,:] = 0
            
            compress.ux[:,:,0] = compress.ux[:,:,-1] = 0
            compress.uy[:,:,0] = compress.uy[:,:,-1] = 0
            compress.uz[:,:,0] = compress.uz[:,:,-1] = 0
            pass

        elif para.boundary_u == 'PFS':
            # Periodic boundary condition for velocity in x and y-direction and free-slip in z-direction
            compress.ux[0,:,:] = compress.ux[-1,:,:] 
            compress.uy[0,:,:] = compress.uy[-1,:,:] 
            compress.uz[0,:,:] = compress.uz[-1,:,:] 

            compress.ux[:,0,:] = compress.ux[:,-1,:] 
            compress.uy[:,0,:] = compress.uy[:,-1,:] 
            compress.uz[:,0,:] = compress.uz[:,-1,:] 
 
            compress.ux[:,:,0] = compress.ux[:,:,1] 
            compress.ux[:,:,-1] = compress.ux[:,:,-2]
            compress.uy[:,:,0] = compress.uy[:,:,1] 
            compress.uy[:,:,-1] = compress.uy[:,:,-2]       
            compress.uz[:,:,0] = compress.uz[:,:,-1] = 0
            pass

        elif para.boundary_u == 'PNS':
            # Periodic boundary condition for velocity in x and y-direction and no-slip in z-direction
            compress.ux[0,:,:] = compress.ux[-1,:,:] 
            compress.uy[0,:,:] = compress.uy[-1,:,:] 
            compress.uz[0,:,:] = compress.uz[-1,:,:] 

            compress.ux[:,0,:] = compress.ux[:,-1,:] 
            compress.uy[:,0,:] = compress.uy[:,-1,:] 
            compress.uz[:,0,:] = compress.uz[:,-1,:] 
            
            compress.ux[:,:,0] = compress.ux[:,:,-1] = 0
            compress.uy[:,:,0] = compress.uy[:,:,-1] = 0
            compress.uz[:,:,0] = compress.uz[:,:,-1] = 0
            pass

        elif para.boundary_u == 'PD':
            # Periodic boundary condition for velocity in x and y-direction and no-slip in z-direction
            compress.ux[0,:,:] = compress.ux[-1,:,:] 
            compress.uy[0,:,:] = compress.uy[-1,:,:] 
            compress.uz[0,:,:] = compress.uz[-1,:,:] 

            compress.ux[:,0,:] = compress.ux[:,-1,:] 
            compress.uy[:,0,:] = compress.uy[:,-1,:] 
            compress.uz[:,0,:] = compress.uz[:,-1,:] 
        
            compress.ux[:,:,0] = compress.ux[:,:,-1] 
            compress.uy[:,:,0] = compress.uy[:,:,-1] 
            compress.uz[:,:,0] = compress.uz[:,:,-1] 
            pass

        pass  

    elif para.dim == 2:

        if para.boundary_u == 'FS':
            # Free-slip boundary condition for velocity at both walls
            compress.ux[0,:] = compress.ux[-1,:] = 0
            compress.uz[0,:] = compress.uz[1,:] 
            compress.uz[-1,:]  = compress.uz[-2,:]

            compress.ux[:,0] = compress.ux[:,1] 
            compress.ux[:,-1] = compress.ux[:,-2]   
            compress.uz[:,0] = compress.uz[:,-1] = 0
            pass
        
        elif para.boundary_u == 'NS':
            # No-slip boundary condition for velocity at both walls
            compress.ux[0,:] = compress.ux[-1,:] = 0
            compress.uz[0,:] = compress.uz[-1,:] = 0
            
            compress.ux[:,0] = compress.ux[:,-1] = 0
            compress.uz[:,0] = compress.uz[:,-1] = 0
            pass

        elif para.boundary_u == 'PFS':
            # Periodic boundary condition for velocity in x-direction and free-slip in z-direction
            compress.ux[0,:] = compress.ux[-1,:]
            compress.uz[0,:] = compress.uz[-1,:]

            compress.ux[:,0] = compress.ux[:,1] 
            compress.ux[:,-1] = compress.ux[:,-2]
            compress.uz[:,0] = compress.uz[:,-1] = 0
            pass

        elif para.boundary_u == 'PNS':
            # Periodic boundary condition for velocity in x-direction and no-slip in z-direction
            compress.ux[0,:] = compress.ux[-1,:]
            compress.uz[0,:] = compress.uz[-1,:]
            
            compress.ux[:,0] = compress.ux[:,-1] = 0
            compress.uz[:,0] = compress.uz[:,-1] = 0
            pass
        
        elif para.boundary_u == 'PD':
            # Periodic boundary condition for velocity in x-direction and z-direction
            compress.ux[0,:] = compress.ux[-1,:]
            compress.uz[0,:] = compress.uz[-1,:]
            
            compress.ux[:,0] = compress.ux[:,-1] 
            compress.uz[:,0] = compress.uz[:,-1] 
            pass

        pass 

    elif para.dim == 1:
        if para.boundary_u == 'NS':
            # No-slip boundary condition for velocity
            compress.uz[0] = compress.uz[-1] = 0
            pass
        elif para.boundary_u == 'PD':
            # Periodic boundary condition for velocity
            compress.uz[0] = compress.uz[-1]
            pass
        pass
        
    pass


def imposeBC_T(compress = Compressible()):
    # Boundary condition for temperature

    if para.dim == 3:
        if para.boundary_T == 'ADF':
            # Fixed boundary on vertical walls, adiabatic on horizontal walls for temperature
            compress.T[0,:,:] = compress.T[1,:,:] 
            compress.T[-1,:,:] = compress.T[-2,:,:]

            compress.T[:,0,:] = compress.T[:,1,:] 
            compress.T[:,-1,:] = compress.T[:,-2,:]

            compress.T[:,:,0] = grid.T_b
            compress.T[:,:,-1] = grid.T_t
            pass

        elif para.boundary_T == 'PDF':
            # Fixed boundary on vertical walls, periodic on horizontal walls for temperature
            compress.T[0,:,:] = compress.T[-1,:,:]

            compress.T[:,0,:] = compress.T[:,-1,:] 

            compress.T[:,:,0] = grid.T_b
            compress.T[:,:,-1] = grid.T_t
            pass

        elif para.boundary_T == 'PD':
            # Periodic boundary conditions on both walls for temperature
            compress.T[0,:,:] = compress.T[-1,:,:]
            compress.T[:,0,:] = compress.T[:,-1,:] 
            compress.T[:,:,0] = compress.T[:,:,-1] 
            pass

        pass

    elif para.dim == 2:
        if para.boundary_T == 'ADF':
            # Fixed boundary on vertical walls, adiabatic on horizontal walls for temperature
            compress.T[0,:] = compress.T[1,:] 
            compress.T[-1,:] = compress.T[-2,:]

            compress.T[:,0] = grid.T_b
            compress.T[:,-1] = grid.T_t
            pass

        elif para.boundary_T == 'PDF':
            # Fixed boundary on vertical walls, periodic on horizontal walls for temperature
            compress.T[0,:] = compress.T[-1,:] 

            compress.T[:,0] = grid.T_b
            compress.T[:,-1] = grid.T_t
            pass

        elif para.boundary_T == 'PD':
            # Periodic boundary conditions on both walls for temperature
            compress.T[0,:] = compress.T[-1,:] 
            compress.T[:,0] = compress.T[:,-1] 
            pass
        pass

    elif para.dim == 1: 
        # Fixed boundary
        compress.T[0] = grid.T_b
        compress.T[-1] = grid.T_t
        pass

    pass

def boundary(compress = Compressible()):
    imposeBC_u(compress)
    imposeBC_T(compress)
    pass