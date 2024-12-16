import os
import pandas
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FFMpegWriter
from scipy.interpolate import griddata

def contour_plot(r_values, theta_values, z_values, var_name, sub_plt=None):
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
    fig, ax = (sub_plt if sub_plt else plt.subplots())
    ax.set_title(var_name)
    contour = ax.contourf(X,Y,Z, levels=50, cmap='viridis')
    if(not sub_plt):
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

def create_movie(r, theta, omega, iterations):
    plt.rcParams['animation.ffmpeg_path'] = "/opt/homebrew/bin/ffmpeg"
    metadata = dict(title='Movie', artist='codinglikemad')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    subplt = fig, ax = plt.subplots()

    fdim = len(r)
    with writer.saving(fig, "omega.mp4", 400):
        def create_frame(it):
            print("creating frame for iteration: " + str(it))
            iteration_index = fdim*(it-1)
            it_omega = omega[iteration_index:iteration_index+fdim]
            contour_plot(r, theta, it_omega, "Omega", subplt)
            writer.grab_frame()
            plt.cla()
        list(map(create_frame, iterations))

def snapshot(r, theta, phi, U_r, U_theta, iteration):
    fdim = len(r)
    iteration_index = fdim*(iteration-1)
    it_omega = omega[iteration_index:iteration_index+fdim]
    it_phi   = phi[iteration_index:iteration_index+fdim]
    it_U_r   = U_r[iteration_index:iteration_index+fdim]
    it_U_theta = U_theta[iteration_index:iteration_index+fdim]

    print("Omega Plot")
    contour_plot(r,theta, it_omega, "Iteration: " + str(iteration) + " Omega")
    print("Phi Plot")
    contour_plot(r,theta, it_phi,   "Iteration: " + str(iteration) + " Phi")
    print("U_r Plot")
    contour_plot_vector(r,theta, it_U_r, it_U_theta, "x", "Iteration: " + str(iteration) + " U_x")
    print("U_theta Plot")
    contour_plot_vector(r,theta, it_U_r, it_U_theta, "y", "Iteration: " + str(iteration) + " U_y")

if __name__ == "__main__":
    iterations = 24

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

    #create_movie(r, theta, omega, range(1,iterations+1))

    #snapshot(r, theta, phi, U_r, U_theta, 1)
    #snapshot(r, theta, phi, U_r, U_theta, 2)
    #snapshot(r, theta, phi, U_r, U_theta, 3)
    snapshot(r, theta, phi, U_r, U_theta, iterations)
    plt.show()
