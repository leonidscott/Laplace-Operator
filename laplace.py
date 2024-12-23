import ctypes
import pathlib
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import csv
from functional import seq
from functools import reduce

def polar_laplace(R_min, R_max, N, bc_type):
    bc_int = (0 if bc_type == "sin" else 1)
    # Setup Interface to C++ Function
    libname=pathlib.Path().absolute() / "numerical-methods.o"
    c_lib = ctypes.CDLL(libname)
    polar_laplace = c_lib.polar_laplace
    polar_laplace.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]
    polar_laplace.restype = ctypes.POINTER(ctypes.c_double)
    # Setup Interface to C++ Function

    arr = polar_laplace(ctypes.c_double(R_min), ctypes.c_double(R_max), N, bc_int)
    print("just finished")
    #for i in range(N*N):
    #    print(arr[i])
    z_values = list(map((lambda i : arr[i]), range(N*N)))

    # Free pointer
    #free = c_lib.freeme
    #free.argtypes = [ctypes.POINTER(ctypes.c_double)]
    #free.restype= None
    #free(arr)

    return z_values


def polar_plot(r_values, t_values, z_values):
    p2c = (lambda r,t : (r*math.cos(t), r*math.sin(t)))
    x_values, y_values = zip(*list(map(p2c, r_values, t_values)))
    #print(x_values)
    #print(y_values)

    print("gets here")
    fig = plt.figure()
    print("but not here")
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(x_values, y_values, z_values, cmap=plt.cm.YlGnBu_r)
    plt.show()

def create_csv(r_points,t_points,z_points):
    print("-- Writing CSV")
    cwd = os.path.abspath(os.getcwd())
    name = cwd + "/laplace-data.csv"
    print("-- Writing CSV to: " + name)
    fields = ['r','t', 'z']
    data = list(map((lambda r,t,z : [r, t, z]), r_points, t_points, z_points))
    with open(name, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(fields)
        writer.writerows(data)
    print("-- CSV Writing Complete")

def d_range(start, end, N, st_ex=False):
    step = (end-start)/(N+1 if st_ex else N)
    np_range = np.arange(
        (start+step if st_ex else start),
        end,
        step
    )
    return np_range[0:N].tolist()

if __name__ == "__main__":
    R_min = 1
    R_max = 5
    N = 6
    bc_type = "sin"
    r_values = seq(d_range(R_min, R_max, N, st_ex=True)) \
        .map(lambda r : [r]*N) \
        .reduce(lambda a,b: a+b) \
        .to_list()
    #r_values += [R_max]*N
    #print(r_values)
    t_values = list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))
    #print(r_values)
    #print(t_values)
    #t_values += t_values[0:N-1]
    #print(t_values)

    print("===== C++ Starts ===== ")
    z_values = polar_laplace(R_min, R_max, N, bc_type)
    print("===== C++ Ends ===== ")
    #print(r_values)
    #print(t_values)
    #bc_pts = list(map((lambda r,t: r * math.sin(t)), r_values[-N:], t_values[-N:]))
    #z_values += bc_pts
    #print(z_values)
    create_csv(r_values, t_values, z_values)

    # Plot Values
    polar_plot(r_values, t_values, z_values)
