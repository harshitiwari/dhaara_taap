import sys
import time 
sys.path.append("lib")
import para as para
from compressible import Compressible
from fns import time_advance_euler, time_advance_RK2, time_advance_RK3
from boundary import boundary
sys.path.append("input")
from init_fields import init_fields

ti = time.time()                    # Initial time

compress = Compressible()           # Making of Compressible object
compress.set_arrays()               # Setting of all arrays
init_fields(compress)               # Initiating initial field profiles at initial time
boundary(compress)                  # Boundary condition for initial field profiles                

if para.Scheme == 'EULER':
    time_advance_euler(compress)    # Time advance steps using Euler method
    pass

elif para.Scheme == 'RK2':
    time_advance_RK2(compress)      # Time advance steps using RK2 method
    pass

elif para.Scheme == 'RK3':
    time_advance_RK3(compress)      # Time advance steps using RK3 method
    pass

tf = time.time()                    # Final time

print('# Total time of simulation = ', tf-ti)        # Time taken to run the code


