from pathlib import Path

main_path = Path(__file__).parents[1].absolute()
path = str(main_path)

type = 'field'                  # Type of plot: field, energy, Nusselt, criteria

Nq = 8                          # Quiver plot distance 
t_start = 0                     # Starting t for field plot

output_folder = '/output/run1/' # Output runs folder
output_file = 'run1.txt'      # Output file name

import field_plots
import output_plots

