import numpy as np
import math
import matplotlib.pyplot as plt

global N  # number of nodes
global dt # timestep
global xl # left boundary
global xr # right boundary
global BC # uxL,uxR,uL,ur 
global ti
global tf

def map_x_to_xi(e, h, xi):
    x_values = [(elem + 1) * (h / 2) + xi for elem in e]
    return x_values

def basis_functions(xi):
    phi1 = 1 - xi / 2
    phi2 = 1 + xi / 2
    return phi1, phi2

def basis_derivatives(i):
    if i == 0:
        return 1/2
    else:
        return -1/2
        
def f(x,t):
    return (pow(math.pi,2) - 1) *  math.exp(-t) * math.sin(math.pi*x) 
def solve_pde(N,dt,xl,xr,BC,ti,tf,ux0, ux1):
    #create uniform grid and connectivity map 
    Ne = N-1
    h = (xr - xl)/Ne
    x = np.linspace(xl,xr,N)
    iee = [[0 for _ in range(2)] for _ in range(Ne)]
    e = [-1/math.sqrt(3), 1/math.sqrt(3)]

    for i in range(Ne):
        x[i]= xl + (i)*h
        iee[i][0] = i
        iee[i][1] = i + 1
    x[N-1]=xr

    z = []
    for x_i in x:
        xi_values = map_x_to_xi(e, h, x_i)
        z.extend([xi_values])

    dxde = h/2
    dedx = 2/h

    K = [[0 for _ in range(N)] for _ in range(N)] # global stiffness
    M = [[0 for _ in range(N)] for _ in range(N)] 
    F = [0] * N # global RHS
    klocal = [[0 for _ in range(2)] for _ in range(2)] # local elements stiffness
    mlocal = [[0 for _ in range(2)] for _ in range(2)]
    flocal = [0] * 2 #local elements RHS

    # Initialize the mass matrix
    M = np.zeros((N, N))
    # Set the main diagonal
    np.fill_diagonal(M, (h * 2)/3)  # 1/6
    # Set the elements above and below the main diagonal
    np.fill_diagonal(M[:, 1:], h/6)  # 1/12
    np.fill_diagonal(M[1:, :], h/6)  # 1/12
    # Set the first and last diagonal elements to 1/6
    M[0, 0] = h / 3
    M[-1, -1] = h / 3
    u = np.sin(np.pi * x)

    for k in range(Ne):# loop over grid elements
        for l in range(2): #local element calculation
            global_node1 = iee[k][l]
            weights = [1,1]
            basis_functions_output = basis_functions(z[k][l])
            phiL,_ = basis_functions(e[0])
            phiR,_ = basis_functions(e[1])
            flocal[l] =  (f(z[k][0],tf) * phiL  + f(z[k][1],tf) * phiR) * dxde

            for m in range(2):
                dphim = basis_derivatives(m)
                dphil = basis_derivatives(l)
                #phi1, phi2 = basis_functions(all_xi[k][l])
                klocal[l][m] = 2*((dphil * dedx) * (dphim * dedx) * dxde)
        #K assembly
        for l in range(2):
            global_node1 = iee[k][l]
            F[global_node1]+= flocal[l]
            for m in range(2):
                global_node2 = iee[k][m]
                K[global_node1][global_node2] += klocal[l][m]
    F[-1] = -F[0]
    
    M_inv = np.linalg.inv(M)
    I = np.identity(N)
    F = np.array(F)
   
    # Dirchilet BC
    K = np.array(K)
    for i in range(N):
        if i == 0 or i == N-1:
            for j in range(N):
                if  j!=i: # set col to 0
                    K[i,:] = 0
                    K[:,i] = 0
                    M[i,:] = 0
                    M[:,i] = 0
                    
                K[0][0] = 1
                K[-1][-1] = 1
                M[0][0] = 1
                M[-1][-1] = 1
                F[i] = K[i][i] * ux0
                F[i] = F[i] - K[i][i] * ux1
    
 

    u = np.dot(np.dot(I - dt * M_inv, K),u) + np.dot(dt * M_inv, F)    
    numerical_sol = u.copy()
    x2 = np.linspace(0,1,100)

    analytical_sol = np.exp(-tf) * np.sin(math.pi * x2)

    plt.plot(x, u, label="Numerical Soln at tf")
    plt.plot(x2, analytical_sol, label="Analytical Soln at tf", linestyle='dashed')
    plt.xlabel("x")
    plt.ylabel("u(x,tf)")
    plt.legend()
    plt.title("Heat Transfer w/ Forward Euler")
    plt.show()

def main():

    N = 11 # number of nodes
    dt = .2 # timestep
    xl = 0 # left boundary
    xr = 1
    BC = [0,0,0,0]#uxL,uxR,uL,ur 
    ti = 0
    tf = 1
    ux0 = 0
    ux1 = 0

    solve_pde(N,dt,xl,xr,BC,ti,tf,ux0,ux1)


if __name__ == "__main__":
    main()
