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

# Forward Euler
Below is the plot generated using forwardEuler.py with number of nodes (N) being 11 and dt=1/551 at the final timestep. As seen in the plot, with this small of a timestep forward euler is stable. 

![plots](plots/goodenoughplot1.png)

I then evaluated the effect of a smaller timestep. While I was unable to determine the exact point at which the solution became unstable, based on trial and error I was able to determine that the solution becomes largely unstable at a timestep of 0.2. The instability is characterized by both ends of the solution not abiding to the boundary conditions.

![plots](plots/unstable_plot0_point_2.png)

I then evaluated the effect that reducing the number of nodes would have on the solution. As seen in the image the accuracy of the solution decreases as we are unable to have enough points to have a smooth curve that resembles the analytical solution.

![plots](plots/4PtForward.png)

# Backward Euler
Below is the plot generated using backwardEuler.py with number of nodes (N) being 11 and dt=1/551 at the final timestep. As seen in the plot, the solution is stable and resembles the stable solution of the forward euler solution. 

![plots](plots/backwrd.png)

I then evaluated the effect of a smaller timestep. While I was unable to determine the exact point at which the solution became unstable, based on trial and error I was able to determine that the solution becomes largely unstable at a timestep of 0.2. The instability is characterized by both ends of the solution not abiding to the boundary conditions.

![plots](plots/unstable_plot0_point_2.png)

I then evaluated the effect that reducing the number of nodes would have on the solution. As seen in the image the accuracy of the solution decreases as we are unable to have enough points to have a smooth curve that resembles the analytical solution.

![plots](plots/4PtForward.png)


