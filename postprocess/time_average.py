import numpy as np
import sys
import output_plots as op

from plot import Path as Path
import plot as pl
sys.path.insert(1, Path)
import para
import grid
from derivative import *

V_rms_av = np.sum(op.V_rms)/(grid.Nf)           # Time average of V_rms

Re = V_rms_av/grid.C2                           # Reynolds number, Re = V_rms_av/(sqrt(Pr/Ra))

print(V_rms_av, Re)