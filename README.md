# Modernized BEM for Time-Domain Elastodynamics (QUADPLET)

This repository contains a modernized Fortran 90 implementation of the **QUADPLET** program, originally developed by José Dominguez in his classic textbook *Boundary Elements in Dynamics*.

The program uses the **Boundary Element Method (BEM)** with quadratic boundary elements to solve 2D (plane) elastic time-domain dynamic problems.

## 🚀 Modernization Features

The original Fortran 77 codebase has been fully refactored to meet modern high-performance computing standards:

1. **Free-Form & Modular Architecture**: The monolithic code was split into clean F90 modules (`data_mod`, `io_mod`, `assembly_mod`, `solver_mod`), eliminating all legacy `COMMON` blocks and enabling strict compiler type checking (`IMPLICIT NONE`).
2. **Dynamic Memory Allocation**: Static limitations (e.g., `PARAMETER(NNE=60)`) were completely removed. The program now reads the exact number of boundary elements and internal points from the input file and allocates arrays dynamically using `ALLOCATABLE`.
3. **OpenMP Parallelization**: The heavy matrix integration routine (`GHMAT12`) that builds the $G$ and $H$ matrices is fully parallelized with OpenMP, mapping independent collocation points to different CPU threads.
4. **LAPACK Solver Integration**: The native manual LU solver (`CSOLVER`) was replaced with LAPACK's highly optimized `ZGESV` (Double-precision Complex General Linear Equation Solver) for massive performance gains in large systems.
5. **Double Precision Upgrade**: The solver inherently enforces **64-bit double precision** using compiler flags (`-fdefault-real-8`), promoting all `REAL` to `REAL*8` and `COMPLEX` to `COMPLEX*16` guaranteeing highly precise numerical integration and results matching out to 18 decimal places.
6. **Python Pre/Post-Processing**: Python scripts (`preproc.py` and `posproc.py`) are provided to integrate with Gmsh, plot geometries, and validate numerical results against analytical solutions.
7. **Colab Ready**: A Jupyter Notebook is included to easily run the entire pipeline in Google Colab.

## 📂 Project Structure

### Fortran Source Code
- `main.f90`: Main program entry point.
- `data_mod.f90`: Global variable declarations and dynamic array allocations.
- `assembly_mod.f90`: BEM influence matrices integration and assembly ($G$ and $H$).
- `io_mod.f90`: Input/Output routines (`INPUT_DATA`, `INTER12`, `OUTPUT12`, `SAVE_DISP`).
- `solver_mod.f90`: System of equations solver wrapper for LAPACK.
- `bessel_mod.f90`: Mathematical routines for computing Bessel functions.

### Python Tools & Validation
- `preproc.py`: Reads Gmsh files (`.msh` via `meshio`), displays the geometry, and exports the simulation input data format (`str12.dat`).
- `posproc.py`: Post-processing script that reads the exported nodal displacements and tractions (`disp_trac.dat`) and plots the numerical BEM solution against the 1D P-wave analytical bar resonance solution.
- `Run_quadplet_on_Colab.ipynb`: Google Colab notebook for running the preprocessing, compiling the BEM code, and executing the analytical validation entirely in the cloud.

### Build Scripts
- `Makefile` / `build.bat`: Build scripts for Windows/MSYS2 and Linux environments.

## ⚙️ Compilation & Build

### Requirements
- **GFortran** (GCC)
- **LAPACK & BLAS** libraries
- **OpenMP** (libgomp)

If you are using Windows with MSYS2, ensure the dependencies are installed:
```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-make
```

### Building the Executable
Navigate to the repository folder and simply run:
```bash
make
```
*(Or `mingw32-make` on MSYS2 Windows)*. This will invoke `gfortran` with `-fopenmp -O3 -fdefault-real-8` and link against `-llapack -lblas` to produce `quadpleh_modern.exe`.

## 🏃 Usage & Validation Pipeline

1. **Preprocessing**: Create your geometry and mesh using Gmsh, and save it as `placa.msh` (or any `.msh`). Run the Python preprocessor to generate the standard input file:
   ```bash
   python preproc.py
   ```
   *Outputs: `str12.dat` and `geometry.png`.*
   
2. **Execution**: Run the Fortran BEM solver:
   ```bash
   ./quadpleh_modern.exe
   ```
   The solver will read `str12.dat`, compute the dynamic matrices across multiple threads, and solve the system via LAPACK. It writes standard outputs to `str12.out` and the specific node tracking data to `disp_trac.dat`.

3. **Post-Processing & Validation**: Compare the BEM numerical results with the analytical solution for infinite strip under uniform step load:
   ```bash
   python posproc.py
   ```
   This script reads the Fortran  outputs  from `disp_trac.dat`, calculates the amplitude $u_y$, compares with the theoretical values, and exports the plot to `result_plot.png`.

## 📚 Reference
- Dominguez, J. (1993). *Boundary Elements in Dynamics*. Computational Mechanics Publications, Southampton.
