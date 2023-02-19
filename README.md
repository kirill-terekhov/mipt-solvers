# mipt-solvers
Supplimentary material for MIPT course "Practical methods for system solutions"

* poisson.cpp     - generates linear system of Poisson problem: -laplace(p) = b on cartesian rectangular NxN grid
* biharmonic.cpp  - generates linear system of Biharmonic problem: -laplace(laplace(p)) = b on cartesian rectangular NxN grid
* stokes.cpp      - generates saddle-point linear system of Stokes problem: -mu laplace(u) + grad(p) = 0, div(u) = 0 on cartesian rectangular NxM grid
* biot.cpp        - generates linear system of Biot problem: -mu laplace(u) - (mu+lambda)*grad(div(u)) + alpha grad(p) = 0, zeta dt(p) + alpha dt(div(u)) - kappa laplace(p) = q.
* deadoil.cpp     - generates system of first Newton step of two-phase filtraion problem: phi dt(S) - div(S^2 kappa grad(p)) = qo, phi dt(1-S) - div((1-S)^2 kappa grad(p)) = qw.
* stokes_grid.cpp - writes solution vector from Stokes or Biot problem to grid.vtk file
* scalar_grid.cpp - writes solution vector from poisson or biharmonic problem to grid.vtk file
* deadoil_grid.cpp - writes combination of initial and solution vectors from deadoil problem to grid.vtk file
