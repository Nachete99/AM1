
# -----------------------------------------------------------------------
#                            IMPORT
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------
#                           FUNCTIONS
# -----------------------------------------------------------------------

def F_Kepler(U_):

    det = (U_[0]**2 + U_[1]**2)**0.5
 
    F = np.array([U_[2], U_[3], -U_[0]/det**3, -U_[1]/det**3 ])

    return F

def Euler(xn, yn, xdn, ydn, dt, n_lim):

    U = np.array(np.zeros((n_lim, 4)))
    U[0] = [xn, yn, xdn, ydn]

    for i in range(len(U[:,1])-1):

        U[i+1] = U[i] + dt*F_Kepler(U[i])
    
    return U[:,0], U[:,1], U[:,2], U[:,3]

def plot(x, y):

    # Create the scatter plot for points only
    plt.plot(x, y,  marker='o', markersize=1, linestyle='-', color='b')  # 's' controls the size of the markers


    # Add labels and a title
    plt.xlabel('x (Position)')
    plt.ylabel('y (Position)')
    plt.title('Vector Field with Starting Points (x, y) and Vector Components (U, V)')

    # Show the plot with a grid
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


# -----------------------------------------------------
#                    MAIN
# -----------------------------------------------------

if __name__ == "__main__":

    xn_0=1; yn_0 = 0; xdn_0 = 0; ydn_0 = 1 # Declaración de condiciones iniciales
    dT = 0.1
    n_lim = 500

    x, y, _, _ = Euler(xn_0, yn_0, xdn_0, ydn_0, dT, n_lim)

    plot(x, y)
