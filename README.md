# DNS-LBM
D3Q19 LBM code for direct numerical simulation (DNS) of turbulence with grid resolution of 133 melion cells. The parallelized code simulates the Taylor-Green vortex at a Reynolds number of 1600, with results matching high-resolution DNS spectral methods.

Solver is written with Intel FORTRAN and parallelized using openMP

# Compile and run

1. Ensure sufficient memory is allocated:
   ```bash
   export OMP_STACKSIZE=6G
2. Compile the code using the provided Makefile
   ```bash
   make
   
3. Run the program: 
   ```bash
   ./scmp_program.

# Computational cost

The LBM simulation presented below, compared with a reference spectral method, was performed using 24 threads on an Intel Xeon E7-8870 processor with a clock speed of 2.4 GHz. The simulation ran for 10 days.

GitPlot.png
