import os
import pandas
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from functional import seq
from functools import reduce
import laplace

def polar3d_plot(r_values, t_values, z_values, title):
    p2c = (lambda r,t : (r*math.cos(t), r*math.sin(t)))
    x_values, y_values = zip(*list(map(p2c, r_values, t_values)))

    fig = plt.figure()
    plt.title(title)
    ax = fig.add_subplot(projection='3d')
    ax.set_ylim(-5,5)
    ax.set_xlim(-5,5)
    ax.plot_trisurf(x_values, y_values, z_values, cmap=plt.cm.YlGnBu_r)

def d_range(start, end, N, st_ex=False):
    step = (end-start)/(N+1 if st_ex else N)
    np_range = np.arange(
        (start+step if st_ex else start),
        end,
        step
    )
    return np_range[0:N].tolist()

def psi_exact(r, t):
    R_max = 5
    R = 1
    U = 1
    # Apply BCs
    if r <= 1:
        return 0
    if r >= R_max:
        return R_max * math.sin(t)
    else:
        return (r - ((pow(R,2)/r))) * U * math.sin(t)

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

def exact_solution(N):
    R_min = 1
    R_max = 5
    r_values = seq(d_range(R_min, R_max, N, st_ex=True)) \
        .map(lambda r : [r]*N) \
        .reduce(lambda a,b: a+b) \
        .to_list()
    t_values = list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))
    z_values = list(map(psi_exact, r_values, t_values))
    return {'r' : r_values, 't' : t_values, 'phi' : z_values}

#Should be in utils
def times(a,b): #Should eventually work with variadic args
    return a*b

if __name__ == "__main__":
    R_min = 1
    R_max = 5
    N = 70
    # Parse Cmd args
    cmd_args = list(map((lambda a : a.lower()), sys.argv))

    # Load Laplace Solver Data
    rs = r_values(R_min, R_max, N)
    ts = t_values(N)
    zs = laplace.polar_laplace(R_min, R_max, N)
    dr = calc_dr(rs)
    dt = ts[1] - ts[0]
    laplace.create_csv(rs, ts, zs)
    csv_name = os.path.abspath(os.getcwd()) + "/laplace-data.csv"
    csv = pandas.read_csv(
        csv_name, skiprows=0,
        dtype={"r" : float, "t" : float, "z" : float}).values.tolist()
    [r_values, t_values, z_values] = list(zip(*csv))
    print("z_values")
    print(max(z_values))
    print(min(z_values))

    z_correction = list(map((lambda z : z*(1.0/(dr*dr))), z_values))
    print("correction")
    print(max(z_correction))
    print(min(z_correction))

    # Exact Solution
    exact = exact_solution(N)
    exact_r, exact_t, exact_phi = exact['r'], exact['t'], exact['phi']
    print("exact")
    print(max(exact_phi))
    print(min(exact_phi))

    print("max error: " + str(max(exact_phi) - max(z_correction)))

    if "Approx3D".lower() in cmd_args:
        polar3d_plot(r_values, t_values, z_correction, "Approx N="+str(N))
    if "Contour".lower() in cmd_args:
        print("Contour in progress")
        #sample = sample_data()
        #r,t,z = sample['r'], sample['t'], sample['z']
        #contour_plot(r_values, t_values, z_values)
        #contour_plot(r, t, z) #TODO: Delete
    if "Exact3D".lower() in cmd_args:
        polar3d_plot(exact_r, exact_t, exact_phi, "Exact N="+str(N))
    if "Error3D".lower() in cmd_args:
        err = list(map((lambda a,b: a-b), exact_phi, z_correction))
        polar3d_plot(r_values, t_values, err, "Error")
    plt.show()
