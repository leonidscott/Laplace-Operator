import os
import pandas
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.patches as patches
from functional import seq
from functools import reduce, partial
import laplace

def polar3d_plot(r_values, t_values, z_values, R_max, title):
    p2c = (lambda r,t : (r*math.cos(t), r*math.sin(t)))
    x_values, y_values = zip(*list(map(p2c, r_values, t_values)))

    fig = plt.figure()
    plt.title(title)
    ax = fig.add_subplot(projection='3d')
    ax.set_ylim(-1 * R_max,R_max)
    ax.set_xlim(-1 * R_max,R_max)
    ax.plot_trisurf(x_values, y_values, z_values, cmap=plt.cm.YlGnBu_r)

def contour_plot(r_values, t_values, z_values, title):
    r = np.array(r_values)
    t = np.array(t_values)
    z = np.array(z_values)

    # Convert to Cartesian Coords
    x = r * np.cos(t)
    y = r * np.sin(t)

    # np grid
    grid_size = len(r_values)
    xi = np.linspace(x.min(), x.max(), grid_size)
    yi = np.linspace(y.min(), y.max(), grid_size)
    X,Y = np.meshgrid(xi, yi)

    # interpolate z values onto grid
    #Z = griddata((x,y), z, (X,Y), method='cubic')

    # Calculate Z Levels
    z_levels = np.linspace(min(z_values), max(z_values), 50)
    z_levels = z_levels[z_levels != 0]

    # Create Plot
    plt.figure()
    plt.title(title)
    #contour = plt.contour(X,Y,Z, levels=z_levels)
    contour = plt.tricontour(x,y,z, levels=z_levels)

    # Add Circle
    circle = patches.Circle((0,0), radius=1, edgecolor='black', facecolor='white')
    plt.gca().add_patch(circle)


def d_range(start, end, N, st_ex=False):
    step = (end-start)/(N+1 if st_ex else N)
    np_range = np.arange(
        (start+step if st_ex else start),
        end,
        step
    )
    return np_range[0:N].tolist()

sin_bc = (lambda R_max, t : R_max * math.sin(t))

def psi_exact(r, t):
    R = 1
    U = 1
    return (r - (pow(R,2)/r)) * U * math.sin(t)

def calc_dr(rs):
    init_val = rs[0]
    for i in range(0,len(rs)):
        if rs[i] != init_val:
            return rs[i] - init_val
    return 0.0

def r_values(R_min, R_max, N):
    return seq(d_range(R_min, R_max, N, st_ex=True)) \
    .map(lambda r : [r]*N) \
    .reduce(lambda a,b: a+b) \
    .to_list()

def t_values(N):
    return list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))

def exact_solution(R_min, R_max, N, bc_type):
    r_values = seq(d_range(R_min, R_max, N, st_ex=True)) \
        .map(lambda r : [r]*N) \
        .reduce(lambda a,b: a+b) \
        .to_list()
    t_values = list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))
    psi = partial(psi_exact)
    z_values = list(map(psi, r_values, t_values))
    return {'r' : r_values, 't' : t_values, 'phi' : z_values}

def call_cpp(R_min, R_max, N, bc_type):
    rs = r_values(R_min, R_max, N)
    ts = t_values(N)
    zs = laplace.polar_laplace(R_min, R_max, N, bc_type)
    dr = calc_dr(rs)
    dt = ts[1] - ts[0]
    laplace.create_csv(rs, ts, zs)

def add_bcs(R_min,R_max,N, r,t,z, bc_type):
    print("len(r): " + str(len(r)))
    print("len(t): " + str(len(t)))
    print("len(z): " + str(len(z)))

    BC_rmin = (lambda t : 0)
    BC_rmax = (psi_exact if bc_type == "exact" else sin_bc)
    # Inner BCs
    inner_r = [R_min]*N
    inner_t = d_range(0, 2*math.pi, N)
    inner_z = list(map(BC_rmin, inner_t))
    # Outer BCs
    outer_r = [R_max]*N
    outer_t = d_range(0, 2*math.pi, N)
    print("len(outer_r): " + str(len(outer_r)))
    print("len(outer_t): " + str(len(outer_t)))
    outer_z = list(map(BC_rmax, outer_r, outer_t))
    # Construct Complete Values
    return {'r' : inner_r + list(r) + outer_r,
            't' : inner_t + list(t) + outer_t,
            'z' : inner_z + list(z) + outer_z}



def print_max_err(exact_vals, approx_vals):
    print("Approximation")
    print("  max: " + str(max(approx_vals)))
    print("  min: " + str(min(approx_vals)))
    print("Exact")
    print("  max: " + str(max(exact_vals)))
    print("  min: " + str(min(exact_vals)))
    print("max error: " + str(max(exact_vals) - max(approx_vals)))

if __name__ == "__main__":
    R_min = 1
    R_max = 15
    N = 20
    bc_type = 'exact'
    # Parse Cmd args
    cmd_args = list(map((lambda a : a.lower()), sys.argv))

    # Load Laplace Solver Data
    if "run".lower() in cmd_args:
        call_cpp(R_min, R_max, N, bc_type)
    csv_name = os.path.abspath(os.getcwd()) + "/laplace-data.csv"
    csv = pandas.read_csv(
        csv_name, skiprows=0,
        dtype={"r" : float, "t" : float, "z" : float}).values.tolist()
    [r_values, t_values, z_values] = list(zip(*csv))

    # Exact Solution
    exact = exact_solution(R_min, R_max, N, bc_type)
    exact_r, exact_t, exact_phi = exact['r'], exact['t'], exact['phi']

    # Add BCs
    aprx_bcs = add_bcs(R_min, R_max, N, r_values, t_values, z_values, bc_type)
    aprx_r, aprx_t, aprx_phi = aprx_bcs['r'], aprx_bcs['t'], aprx_bcs['z']
    print("len(aprx_r): "   + str(len(aprx_r)))
    print("len(aprx_t): "   + str(len(aprx_t)))
    print("len(aprx_phi): " + str(len(aprx_phi)))

    exact_bcs = add_bcs(R_min, R_max, N, exact_r, exact_t, exact_phi, bc_type)
    exact_r, exact_t, exact_phi = exact_bcs['r'], exact_bcs['t'], exact_bcs['z']


    print_max_err(exact_phi, z_values)
    if "AllPlots".lower() in cmd_args:
        cmd_args = ["approx3d", "approxcontour", "exact3d", "exactcontour", "error3d"]
    if "Approx3D".lower() in cmd_args:
        polar3d_plot(aprx_r, aprx_t, aprx_phi, R_max, "Approx N="+str(N))
    if "ApproxContour".lower() in cmd_args:
        contour_plot(aprx_r, aprx_t, aprx_phi, "Approx N="+str(N))
    if "Exact3D".lower() in cmd_args:
        polar3d_plot(exact_r, exact_t, exact_phi, R_max, "Exact N="+str(N))
    if "ExactContour".lower() in cmd_args:
        contour_plot(exact_r, exact_t, exact_phi, "Exact N="+str(N))
    if "Error3D".lower() in cmd_args:
        err = list(map((lambda a,b: a-b), exact_phi, aprx_phi))
        polar3d_plot(aprx_r, aprx_t, err, R_max, "Error")
    plt.show()
