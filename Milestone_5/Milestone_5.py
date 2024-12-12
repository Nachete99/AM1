# -----------------------------------------------------------------------
#                            IMPORT
# -----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import    norm

# -----------------------------------------------------------------------
#                           FUNCTIONS
# -----------------------------------------------------------------------

def N_Body_Problem(U, Nb, Nc):

    pU = np.reshape(U, (Nb, Nc, 2))

    r = np.reshape(pU[:, :, 0], (Nb, Nc))
    v = np.reshape(pU[:, :, 1], (Nb, Nc))

    Fs = np.zeros(2*Nb*Nc)
    pF = np.reshape(Fs, (Nb, Nc, 2))

    drdt = np.reshape(pF[:, :, 0], (Nb, Nc))
    dvdt = np.reshape(pF[:, :, 1], (Nb, Nc))


    for n1 in range(Nb):
        drdt[n1,:] = v[n1,:]
        for n2 in range(Nb):
            if n1 != n2:
                d = r[n2,:] - r[n1,:]
                dvdt[n1,:] = dvdt[n1,:] + d/norm(d)**3

    return Fs


def F_N_Body_Problem(U, t):
    return N_Body_Problem(U, Nb, Nc)


def Cauchy(U_, dT, scheme, F):
    
    for i in range(len(U_[:,1])-1):

        T = dT*(i+1)
        U_[i+1] = scheme(U_[i], T, dT, F)

    return U_


def RK4(U_n, t, dt, F):
    # Application of the Runge Kutta of 4th order numerical scheme
    # Definition of the parameters for Runge Kutta 4rd order

    K1 = F(U_n, t)
    K2 = F(U_n + 0.5*dt*K1, t + 0.5*dt*K1)
    K3 = F(U_n + 0.5*dt*K2, t + 0.5*dt*K2)
    K4 = F(U_n + dt*K3, t + dt)

    return U_n + dt*(K1 + 2*K2 + 2*K3 + K4)/6


# -----------------------------------------------------------------------
#                    MAIN
# -----------------------------------------------------------------------

if __name__ == "__main__":

    dT = 0.001
    N = 3000 # Numero de nodos
    
    # Initial conditions:

    r0=[];v0=[]

    # Cuerpo 1
    r0.append([.5, 0, 0])
    v0.append([0, 1, 0])

    # Cuerpo 2
    r0.append([ -.5, 0, 0])
    v0.append([ 0, -1, 0])

    # Cuerpo 3
    r0.append([ 0, 0, .5 ])
    v0.append([ -1, 0., 0. ])

    # Cuerpo 4
    r0.append([ 0, 0, -.5 ]) 
    v0.append([ 1, 0., 0. ])

    Nb = len(r0)
    Nc = len(r0[0])

    U_ = np.zeros((N, 2*Nc*Nb))
    U0 = np.reshape(U_[0], (Nb, Nc, 2))
    r0_ = np.reshape(U0[:, :, 0], (Nb, Nc))
    v0_ = np.reshape(U0[:, :, 1], (Nb, Nc))
    r0_[:Nb,:] = r0
    v0_[:Nb,:] = v0

    U = Cauchy(U_, dT, RK4, F_N_Body_Problem)

    # Representacion de resultados
    Us  = np.reshape( U, (N, Nb, Nc, 2) ) # Us[tiempo][cuerpo][coordenada][0, 1]
    r   = np.reshape( Us[:, :, :, 0], (N, Nb, Nc) ) 

    
# Representacion de resultados (en 2 dimensiones)
    for i in range(Nb):
        plt.plot(  r[:, i, 0], r[:, i, 1] )
    plt.axis('equal')
    plt.grid()
    plt.show()

    # Representacion de resultados (en 3 dimensiones)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(Nb):
        ax.plot3D(r[:, i, 0], r[:, i, 1], r[:, i, 2], label=f'Cuerpo {i+1}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()