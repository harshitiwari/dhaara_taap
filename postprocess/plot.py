from pathlib import Path

main_path = Path(__file__).parents[1].absolute()
path = str(main_path)

type = 'Nusselt'          # Type of plot: field, field_diff, horiz_profile, 
                        # energy, Nusselt, criteria, ratio, V_rms
                        # T_ts, ux_ts, uz_ts

Nq = 8                          # Quiver plot distance 
scale = 0.5                     # Scale for quiver arrows
t_start = 0                     # Starting t for field plot
point = 'iii'                     # Time series at point i, ii, iii, iv

sec = '1'                       # Section of run
run = '1'                       # Run number

output_folder = '/output/sec' + sec + '/run' + run + '/'       # Output runs folder
output_file = 'sec' + sec + 'run' + run + '.txt'               # Output file name

import field_plots
import output_plots
import timeseries_plots

