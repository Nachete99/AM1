# -----------------------------------------------------------------------
#                            IMPORT
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

# -----------------------------------------------------------------------
#                           FUNCTIONS
# -----------------------------------------------------------------------

def Kepler(U_, t):
    # print(U_)
    det = (U_[0]**2 + U_[1]**2)**0.5
 
    F = np.array([U_[2], U_[3], -U_[0]/det**3, -U_[1]/det**3 ])

    return F


def RK4(U_n, t, dt, F):
    # Application of the Runge Kutta of 4th order numerical scheme
    # Definition of the parameters for Runge Kutta 4rd order

    K1 = F(U_n, t)
    K2 = F(U_n + 0.5*dt*K1, t + 0.5*dt*K1)
    K3 = F(U_n + 0.5*dt*K2, t + 0.5*dt*K2)
    K4 = F(U_n + dt*K3, t + dt)

    return U_n + dt*(K1 + 2*K2 + 2*K3 + K4)/6


def Euler(U_n, t, dt, F):
    # Application of the Euler numerical scheme

        return U_n + dt*F(U_n, dt)


def Crank_Nicolson(U_n, t, dt, F):
    # Application of the Crank Nicolson numerical scheme

    def G(X):
        return X - U_n - dt/2*(F(U_n,t) + F(X,t))

    return newton(G, U_n)


def Inverse_Euler(U_n, t, dt, F):
    # Application of the inverse Euler numerical scheme

    def G(X):
        return X - U_n - dt*F(X,t)
    
    return newton(G, U_n, maxiter=10000)

def Richardson(U_n, t, dt, F):

    return


def Cauchy(U_, dT, scheme, F):
    
    for i in range(len(U_[:,1])-1):

        T = dT*(i+1)
        U_[i+1] = scheme(U_[i], T, dT, F)

    return U_


def Cauchy_Error(U_, dT, scheme, F, output_flag = False, plot_flag = False):

    N = len(U_[:,1])
    U_refined = np.array(np.zeros((2*N-1, 4)))
    U_refined[0, :] = U_[0,:].copy()

    
    U1 = Cauchy(U_, dT, scheme, F)
    U2 = Cauchy(U_refined, dT/2, scheme, F)
    
    Error = np.array(np.zeros((N, 4)))
    t = np.array(np.zeros(N))

    for i in range(0, N):
        t[i] = dT*i
        Error[i,:] = U2[2*i,:] - U1[i,:]

    if plot_flag:
        plot(U1[:,0], U1[:,1], scheme.__name__)
        plot( t, Error[:,0],"Error " + scheme.__name__)
    
    if output_flag:
        return Error
    
    return


def Convergence(U_0, dT, scheme, F):

    Np = 15 #Number of regression points
    logE = np.zeros(Np)
    logN = np.zeros(Np)
    U_ = U_0.copy()
    CI = U_0[0,:]
    N = len(U_[:,1])
    dt = dT

    for i in range(Np):
        
        E = Cauchy_Error(U_, dt, scheme, F, True, False)
        logE[i] = np.log10(np.linalg.norm(E[-1,:]))
        logN[i] = np.log10(N)
        N = 2*N-1
        U_ = np.array(np.zeros((N, 4)))
        U_[0,:] = CI
        dt = 0.5*dt

    plot(logN, logE, f"Convergence graph for {scheme.__name__} scheme")

def plot(x, y, method):

    # Define points (x, y) where vectors start
    x = np.array(x)
    y = np.array(y)

    # Create the scatter plot for points only
    # plt.scatter(x, y, color='b', marker='o', s=1)  # 's' controls the size of the markers
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(x, y, marker='o', linestyle='-', markersize=2, color='b')

    # Add labels and a title
    plt.xlabel('x (Position)')
    plt.ylabel('y (Position)')
    plt.title(f'Kepler solutions to {method} method')

    # Show the plot with a grid
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


# -----------------------------------------------------------------------
#                    MAIN
# -----------------------------------------------------------------------

if __name__ == "__main__":

    xn_0=1; yn_0 = 0; xdn_0 = 0; ydn_0 = 1 # Declaración de condiciones iniciales
    dT = 0.001
    N = 3000 # Numero de nodos

    U = np.array(np.zeros((N, 4)))
    U[0] = [xn_0, yn_0, xdn_0, ydn_0]


    # Cauchy_Error(U, dT, Euler, Kepler)
    # Cauchy_Error(U, dT, RK4, Kepler)
    # Cauchy_Error(U, dT, Crank_Nicolson, Kepler)
    # Cauchy_Error(U, dT, Inverse_Euler, Kepler)

    Convergence(U, dT, Euler, Kepler)