from grid import ncp,copy
import grid

import sys
sys.path.append("..")

import para

# Central-difference method

def dfx_c(f):
    # Derivative of array 'f' w.r.t. x

    # For peridoic f
    if para.boundary_u == 'PFS' or 'PNS' and para.boundary_T == 'PD':
        if para.order == 2:
            # 2nd order accuracy
            grid.temp_d[1:-1,:] = (f[2:,:] - f[:-2,:])/(2*grid.dx)
            # Boundary terms when f is periodic in x
            grid.temp_d[0,:] = (f[1,:] - f[-2,:])/(2*grid.dx)
            grid.temp_d[-1,:] = copy(grid.temp_d[0,:])
            pass

        elif para.order == 4:
            # 4th order accuracy
            grid.temp_d[2:-2,:] = (-f[4:,:] + 8*f[3:-1,:] - 8*f[1:-3,:] + f[:-4,:])/(12*grid.dx)
            # Boundary terms when f is periodic in x
            grid.temp_d[0,:] = (-f[2,:] + 8*f[1,:] - 8*f[-2,:] + f[-3,:])/(12*grid.dx)
            grid.temp_d[1,:] = (-f[3,:] + 8*f[2,:] - 8*f[0,:] + f[-2,:])/(12*grid.dx)
            grid.temp_d[-1,:] = copy(grid.temp_d[0,:])
            grid.temp_d[-2,:] = (f[-4,:] - 8*f[-3,:] + 8*f[0,:] - f[1,:])/(12*grid.dx)
            pass
        pass
    
    # For non-peridoic f
    else:
        if para.order == 2:
            # 2nd order accuracy
            grid.temp_d[1:-1,:] = (f[2:,:] - f[:-2,:])/(2*grid.dx)
            # Boundary terms
            grid.temp_d[0,:] = (-f[2,:] + 4*f[1,:] - 3*f[0,:])/(2*grid.dx)
            grid.temp_d[-1,:] = (f[-3,:] - 4*f[-2,:] + 3*f[-1,:])/(2*grid.dx)
            pass

        elif para.order == 4:
            # 4th order accuracy
            grid.temp_d[2:-2,:] = (-f[4:,:] + 8*f[3:-1,:] - 8*f[1:-3,:] + f[:-4,:])/(12*grid.dx)
            # Boundary terms
            grid.temp_d[0,:] = (-3*f[4,:] + 16*f[3,:] - 36*f[2,:] + 48*f[1,:] - 25*f[0,:])/(12*grid.dx)
            grid.temp_d[1,:] = (-3*f[5,:] + 16*f[4,:] - 36*f[3,:] + 48*f[2,:] - 25*f[1,:])/(12*grid.dx)
            grid.temp_d[-1,:] = (3*f[-5,:] - 16*f[-4,:] + 36*f[-3,:] - 48*f[-2,:] + 25*f[-1,:])/(12*grid.dx)
            grid.temp_d[-2,:] = (3*f[-6,:] - 16*f[-5,:] + 36*f[-4,:] - 48*f[-3,:] + 25*f[-2,:])/(12*grid.dx)
            pass
        pass

    return grid.temp_d

def dfz_c(f):
    # Derivative of array 'f' w.r.t. z

    if para.dim == 2:
        if para.order == 2:
            # 2nd order accuracy
            grid.temp_d[:,1:-1] = (f[:,2:] - f[:,:-2])/(2*grid.dz)
            # Boundary terms
            grid.temp_d[:,0] = (4*f[:,1] - f[:,2] - 3*f[:,0])/(2*grid.dz)
            grid.temp_d[:,-1] = (3*f[:,-1] - 4*f[:,-2] + f[:,-3])/(2*grid.dz)
            pass

        elif para.order == 4:
            # 4th order accuracy
            grid.temp_d[:,2:-2] = (-f[:,4:] + 8*f[:,3:-1] - 8*f[:,1:-3] + f[:,:-4])/(12*grid.dz)
            # Boundary terms
            grid.temp_d[:,0] = (-3*f[:,4] + 16*f[:,3] - 36*f[:,2] + 48*f[:,1] - 25*f[:,0])/(12*grid.dz)
            grid.temp_d[:,1] = (-3*f[:,5] + 16*f[:,4] - 36*f[:,3] + 48*f[:,2] - 25*f[:,1])/(12*grid.dz)
            grid.temp_d[:,-1] = (3*f[:,-5] - 16*f[:,-4] + 36*f[:,-3] - 48*f[:,-2] + 25*f[:,-1])/(12*grid.dz)
            grid.temp_d[:,-2] = (3*f[:,-6] - 16*f[:,-5] + 36*f[:,-4] - 48*f[:,-3] + 25*f[:,-2])/(12*grid.dz)
            pass
        pass

    elif para.dim == 1:
        if para.order == 2:
            # 2nd order accuracy
            grid.temp_d[1:-1] = (f[2:] - f[:-2])/(2*grid.dz)
            # Boundary terms
            grid.temp_d[0] = (4*f[1] - f[2] - 3*f[0])/(2*grid.dz)
            grid.temp_d[-1] = (3*f[-1] - 4*f[-2] + f[-3])/(2*grid.dz)
            pass

        elif para.order == 4:
            # 4th order accuracy
            grid.temp_d[2:-2] = (-f[4:] + 8*f[3:-1] - 8*f[1:-3] + f[:-4])/(12*grid.dz)
            # Boundary terms
            grid.temp_d[0] = (-3*f[4] + 16*f[3] - 36*f[2] + 48*f[1] - 25*f[0])/(12*grid.dz)
            grid.temp_d[1] = (-3*f[5] + 16*f[4] - 36*f[3] + 48*f[2] - 25*f[1])/(12*grid.dz)
            grid.temp_d[-1] = (3*f[-5] - 16*f[-4] + 36*f[-3] - 48*f[-2] + 25*f[-1])/(12*grid.dz)
            grid.temp_d[-2] = (3*f[-6] - 16*f[-5] + 36*f[-4] - 48*f[-3] + 25*f[-2])/(12*grid.dz)
            pass
        pass

    return grid.temp_d
