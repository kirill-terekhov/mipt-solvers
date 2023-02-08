# mipt-solvers
Supplimentary material for MIPT course "Practical methods for system solutions"

* poisson.cpp     - generates linear system of Poisson problem: -laplace(p) = b on cartesian rectangular NxN grid
* biharmonic.cpp  - generates linear system of Biharmonic problem: -laplace(laplace(p)) = b on cartesian rectangular NxN grid
* stokes.cpp      - generates saddle-point linear system of Stokes problem: -mu laplace(u) + grad(p) = 0, div(u) = 0 on cartesian rectangular NxM grid
* biot.cpp        - generates linear system of Biot problem: -mu laplace(u) - (mu+lambda)*grad(div(u)) + alpha grad(p) = 0, zeta dt(p) + alpha dt(div(u)) - kappa laplace(p) = q.
* stokes_grid.cpp - writes solution vector from Stokes or Biot problem to grid.vtk file
* scalar_grid.cpp - writes solution vector from poisson or biharmonic problem to grid.vtk file
