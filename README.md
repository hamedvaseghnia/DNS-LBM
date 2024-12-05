# DNS-LBM
D3Q19 LBM code for direct numerical simulation (DNS) of turbulence. The parallelized code simulates the Taylor-Green vortex at a Reynolds number of 1600, with results matching high-resolution DNS spectral methods



# Compile and run

Ensure sufficient memory is allocated: export OMP_STACKSIZE=6G
Compile the code using the provided Makefile: make
Run the program: ./scmp_program.
