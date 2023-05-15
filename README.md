# dhaara_taap
## Compressible Convection Solver

- This is a sequential Object Oriented Python code which solves fully compressible convection in Cartesian box. 
- This is GPU enabled. 
- Currently, I have implemented it for 1D and 2D.

Steps to run the code:
1. Go to `para.py` file for changing the parameters. 
  - To make it run on GPU, put GPU for `device` and its rank. 
  - Make a output directory where you want to store the fields and give its path to `output_dir`. 
  - You can change the time advance scheme from `Scheme`.
  - Put the grid parameters and control parameters.
2. To start the solver run `main.py` using Python. You can save the output in a text file using one of the following commands:
  - `python3 main.py > output\output.txt`
  - `nohup python3 main.py &> output\output.txt &`
3. You can postprocess the output from by running the `postprocess\plot.py` file.
  - Put the output fields folder directory in `output_folder` and output file name in `output_file`.
  - Change the `type` and run `plot.py` file to obtain various graphs.

About the output fields and output print parameters:
- The output field contains: rho, ux, uz, T and parameters. 
- T and rho are actually the perturbations from the equilibrium profiles, T-T_0 and rho - rho_0.
- parameters contain: alpha, m, epsilon, Pr, Ra, dt in sequence.
- The output print contains all the parameters of the files in first three lines with comment #.
- The next lines contain the following columns in order: t, Ke, Ie, V_rms, Nu_b, Nu_t, cri, p_t, g_t.
- Ke = Total kinetic energy, Ie = Total internal energy, V_rms = Root mean square of velocity volume average, Nu_b, Nu_t = Nusselt number at bottoma and top plate, cri = difference in Schwarzschild criteria, p_t, g_t = Volume average of pressure gradient term and gravity term in uz equation respectively.
 

Notes:
- The solver needs the initial fields. I will recommend to leave them since they contain the equilribrium fields. 
- If you want to change them, you can give them from `lib\compressible.py` file with the function `init_fields`. 
