
# -----------------------------------------------------------------------
#                            IMPORT
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------
#                           FUNCTIONS
# -----------------------------------------------------------------------

def F(x, y, x_dot, y_dot):

    F1 = x_dot
    F2 = y_dot
    F3 = -x/(x**2 + y**2)**(3/2)
    F4 = -y/(x**2 + y**2)**(3/2) 

    return F1, F2, F3, F4

def U(xn, yn, xdn, ydn, n_lim):

    x = [xn]; y = [yn]; U = [xdn]; V = [ydn]

    n = 0
    while n <= n_lim:

        F1, F2, F3, F4 = F(x[n], y[n], U[n], V[n])

        x.append(x[n] + dT * F1)
        y.append(y[n] + dT * F2)
        U.append(U[n] + dT * F3)
        V.append(V[n] + dT* F4)
        
        n += 1

    return x, y, U, V

def plot(x, y, U, V, vector_flag):

    # Define points (x, y) where vectors start
    x = np.array(x)
    y = np.array(y)

    # Define vector components (U, V) for each point
    U = np.array(U)  # X components of vectors
    V = np.array(V)  # Y components of vectors

    if vector_flag:
        # Create the quiver plot
        plt.quiver(x, y, U, V, angles='xy', scale_units='xy', scale=1, color='b')
    else:
        # Create the scatter plot for points only
        plt.scatter(x, y, color='b', marker='o', s=1)  # 's' controls the size of the markers



    # Set limits for x and y axes
    # plt.xlim(min(x)+2, max(x)+2)
    # plt.ylim(min(y)+2, max(y)+2)

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

xn=1; yn = 0; xdn = 0; ydn = 1
dT = 0.1
n_lim = 50000
plot_vector_flag = False

x, y, U, V = U(xn, yn, xdn, ydn, n_lim)

plot(x, y, U ,V, plot_vector_flag)
