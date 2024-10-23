# -----------------------------------------------------------------------
#                            IMPORT
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------
#                           FUNCTIONS
# -----------------------------------------------------------------------



def Euler(xn, yn, xdn, ydn, dt, n_lim):

    x = [xn]; y = [yn]; dx = [xdn]; dy = [ydn]

    n = 0

    while n <= n_lim:

        F1 = dx[n]
        F2 = dy[n]
        F3 = -x[n]/(x[n]**2 + y[n]**2)**(3/2)
        F4 = -y[n]/(x[n]**2 + y[n]**2)**(3/2)

        x.append(x[n] + dt * F1)
        y.append(y[n] + dt * F2)
        dx.append(dx[n] + dt * F3)
        dy.append(dy[n] + dt * F4)
        
        n += 1
    
    plot(x, y, 'Euler')

    return

def Crank_Nicolson(xn, yn, xdn, ydn, dt, n_lim):

    #TBD

    return

def Runge_Kutta_4(xn, yn, xdn, ydn, dt, n_lim):

    x = [xn]; y = [yn]; dx = [xdn]; dy = [ydn]

    n = 0

    while n <= n_lim:
        
        K1_x = dx[n]
        K1_y = dy[n]
        r1 = (x[n]**2 + y[n]**2)
        K1_dx = -x[n]/r1**(3/2)
        K1_dy = -y[n]/r1**(3/2)
        
        K2_x = dx[n] + 0.5*dt*K1_dx
        K2_y = dy[n] + 0.5*dt*K1_dy
        r2 = ((x[n] + 0.5*dt*K1_x)**2 + (y[n] + 0.5*dt*K1_y)**2)
        K2_dx = -(x[n] + 0.5*dt*K1_x)/r2**(3/2)
        K2_dy = -(y[n] + 0.5*dt*K1_y)/r2**(3/2)

        K3_x = dx[n] + 0.5*dt*K2_dx
        K3_y = dy[n] + 0.5*dt*K2_dy
        r3 = ((x[n] + 0.5*dt*K2_x)**2 + (y[n] + 0.5*dt*K2_y)**2)
        K3_dx = -(x[n] + 0.5*dt*K2_x)/r3**(3/2)
        K3_dy = -(y[n] + 0.5*dt*K2_y)/r3**(3/2)

        K4_x = dx[n] + dt*K3_dx
        K4_y = dy[n] + dt*K3_dy
        r4 = ((x[n] + dt*K3_x)**2 + (y[n] + dt*K3_y)**2)
        K4_dx = -(x[n] + dt*K3_x)/r4**(3/2)
        K4_dy = -(y[n] + dt*K3_y)/r4**(3/2)

        x.append(x[n] + dt * (K1_x + 2*K2_x + 2*K3_x + K4_x)/6 )
        y.append(y[n] + dt * (K1_y + 2*K2_y + 2*K3_y + K4_y)/6 )
        dx.append(dx[n] + dt * (K1_dx + 2*K2_dx + 2*K3_dx + K4_dx)/6 )
        dy.append(dy[n] + dt * (K1_dy + 2*K2_dy + 2*K3_dy + K4_dy)/6 )

        n+=1

    plot(x, y, 'RK4')
        
    return


def Euler_inverso(U, dt, t, F):
    def G(X):
        return  X - U - dt * F(X,t)

    return 


def plot(x, y, method):

    # Define points (x, y) where vectors start
    x = np.array(x)
    y = np.array(y)

    # Create the scatter plot for points only
    plt.scatter(x, y, color='b', marker='o', s=1)  # 's' controls the size of the markers

    # Add labels and a title
    plt.xlabel('x (Position)')
    plt.ylabel('y (Position)')
    plt.title(f'Kepler solutions to {method} method')

    # Show the plot with a grid
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


# -----------------------------------------------------
#                    MAIN
# -----------------------------------------------------

if __name__ == "__main__":
    x_0=1; y_0 = 0; dx_0 = 0; dy_0 = 1 # Condiciones iniciales
    dT = 0.1
    n_lim = 500
    plot_vector_flag = True

    Euler(x_0, y_0, dx_0, dy_0, dT, n_lim)
    Crank_Nicolson(x_0, y_0, dx_0, dy_0, dT, n_lim)
    Runge_Kutta_4(x_0, y_0, dx_0, dy_0, dT, n_lim)
    Euler_inverso()

    print('END OF SCRIPT')
