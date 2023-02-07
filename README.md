# mipt-solvers
Supplimentary material for MIPT course "Practical methods for system solutions"

* poisson.cpp     - generates linear system of -\Delta p = b on cartesian rectangular NxN grid
* biharmonic.cpp  - generates linear system of -\Delta^2 p = b on cartesian rectangular NxN grid
* stokes.cpp      - generates saddle-point linear system of Stokes problem: -mu \Delta u + \nabla p = 0, \nabla^T u = 0 on cartesian rectangular NxM grid
* stokes_grid.cpp - writes solution vector from stokes problem to grid.vtk file
