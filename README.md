# Numerical calculation of velocity profile of a steady flow through a pipe with hexagonal cross-cection (FORTRAN)

## Overview
This project numerically solves the Poisson equation to compute the velocity profile of a Newtonian fluid flowing through a hexagonal pipe cross-section under steady-state laminar flow conditions.

The solution uses a finite-difference method on a custom hexagonal mesh.

## Project History
This project is part of final assignment of course *Mathematical methods for physical chemistry* (Autumn 2015) from MSc studies at the UCT Prague.

## Project Structure
- `src/HEXPF.for`: Main Fortran 77 code implementing the numerical solver.
- `data/velocity_profile.dat`: Computed velocity profile data.
- `data/profile.plt`: Plotting script for visualization with Gnuplot.
- `report/projekt.pdf`: Original project description in Czech (background material).

## How to Run
1. Compile the Fortran code using a compiler such as `gfortran`.
    ```bash
    gfortran -o pipe_flow src/HEXPF.for
    ./pipe_flow
    ```
2. The output `velocity_profile.dat` can be visualized with the provided `profile.plt` script.:
    ```bash
    gnuplot data/profile.plt
    ```

## Requirements
- FORTRAN 77 compatible compiler (e.g., `gfortran`)
- gnuplot (optional, for visualization)

## Notes
- The model description is in Czech, but the mathematical model is the standard finite-difference method.
- The project demonstrates numerical modeling, relaxation methods, and custom mesh generation.

## License
For educational purposes only.