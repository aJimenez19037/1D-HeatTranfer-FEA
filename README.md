# 1D-HeatTranfer-FEA
The purpose of this project was to solve a 1D heat transfer problem given initial and Dirichlet boundary conditions:
```
#Problem
U_t - U_xx = F(x,t) (x,t): (0,1) x (0,1)
F(x,t) = (pi**2 - 1) * e^(-t) * sin(pi * x)
#Dirchelet boundary conditions
U(x,0) = sin(pi * x)
U(0,t) = U(1,t) = 0
```
After solving the problem numerically it was plotted against the numerical solution:
```
U(x,t) = e^(-t) sin(pi * x)
```
# Included Files
## forwardEuler.py
This file contains the solution using forward euler.
## backwardEuler.py 
This file contains the solution using backward euler. Main difference from forward Euler.py is how U is calculated.
## plot
This directory is where different plot are stored. Plots contain different number of nodes and timestep size
## hand-written-work.pdf
This pdf contains my handwritten work which includes weak form derivation, forward euler implementation, gelerkin expansion implementation, calculating mass and stiffness matrix. 



