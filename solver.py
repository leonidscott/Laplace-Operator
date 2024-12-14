import os
import pandas
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import griddata

def contour_plot(r_values, theta_values, z_values, var_name):
    r     = np.array(r_values)
    theta = np.array(theta_values)
    z     = np.array(z_values)

    # Convert to Cartesian Coords
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # np grid
    grid_size = len(r_values)
    xi = np.linspace(x.min(), x.max(), grid_size)
    yi = np.linspace(y.min(), y.max(), grid_size)
    X,Y = np.meshgrid(xi, yi)

    # interpolate z values onto grid
    Z = griddata((x,y), z, (X,Y))

    # Create Plot
    fig, ax = plt.subplots()
    ax.set_title(var_name)
    contour = ax.contourf(X,Y,Z, levels=50, cmap='viridis')
    fig.colorbar(contour, label=var_name)

    # Add Circle
    circle = patches.Circle((0,0), radius=1, edgecolor='black', facecolor='white')
    ax.add_patch(circle)

# For Radial Z Values
def contour_plot_vector(r_values, theta_values, zr_values, ztheta_values, direction, var_name):
    r     = np.array(r_values)
    theta = np.array(theta_values)
    zr    = np.array(zr_values)
    zt    = np.array(ztheta_values)

    # Convert to Cartesian Coords
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    #zx = zr * np.cos(zt)
    #zy = zr * np.sin(zt)
    zx = zr* np.cos(theta) - (r* zt * np.sin(theta))
    zy = zr* np.sin(theta) + (r* zt * np.cos(theta))

    # np grid
    grid_size = len(r_values)
    xi = np.linspace(x.min(), x.max(), grid_size)
    yi = np.linspace(y.min(), y.max(), grid_size)
    X,Y = np.meshgrid(xi, yi)

    # interpolate z values onto grid
    z_dir = (zx if direction == "x" else zy)
    Z = griddata((x,y), z_dir, (X,Y))

    # Create Plot
    fig, ax = plt.subplots()
    ax.set_title(var_name)
    contour = ax.contourf(X,Y,Z, levels=50, cmap='rainbow')
    fig.colorbar(contour, label=var_name)

    # Add Circle
    circle = patches.Circle((0,0), radius=1, edgecolor='black', facecolor='white')
    ax.add_patch(circle)

if __name__ == "__main__":
    iteration = 1

    grid_csv_name = os.path.abspath(os.getcwd()) + "/grid.csv"
    grid_csv = pandas.read_csv(
        grid_csv_name, skiprows=0,
        dtype={"r" : float, "theta" : float}).values.tolist()
    [r, theta] = list(zip(*grid_csv))
    solver_csv_name = os.path.abspath(os.getcwd()) + "/solver-data.csv"
    solver_csv = pandas.read_csv(
        solver_csv_name, skiprows=0,
        dtype={"r" : float, "theta" : float, "omega" : float,
               "phi" : float, "U_r" : float, "U_theta" : float}).values.tolist()
    [omega, phi, U_r, U_theta] = list(zip(*solver_csv))

    fdim = len(r)
    it_omega = omega[iteration-1:iteration-1+fdim]
    it_phi   = phi[iteration-1:iteration-1+fdim]
    it_U_r   = U_r[iteration-1:iteration-1+fdim]
    it_U_theta = U_theta[iteration-1:iteration-1+fdim]

    contour_plot(r,theta, it_omega, "Omega")
    contour_plot(r,theta, it_phi,   "Phi")
    contour_plot_vector(r,theta, it_U_r, it_U_theta, "x", "U_x")
    contour_plot_vector(r,theta, it_U_r, it_U_theta, "y", "U_y")
    plt.show()
