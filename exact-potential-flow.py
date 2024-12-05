import math
import numpy as np
from functional import seq
from functools import reduce
import matplotlib.pyplot as plt

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
        return (r - (R*R/r)) * U * math.sin(t)

def polar_plot(r_values, t_values, z_values):
    p2c = (lambda r,t : (r*math.cos(t), r*math.sin(t)))
    x_values, y_values = zip(*list(map(p2c, r_values, t_values)))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(x_values, y_values, z_values, cmap=plt.cm.YlGnBu_r)
    plt.show()


def contour_plot(r_values, t_values, z_values):
    npr = np.array(list(set(r_values)))
    print(npr)
    npt = np.array(list(set(t_values)))
    print("len(npr): " + str(len(npr)))
    print("len(npt): " + str(len(npt)))
    npz = np.array([z_values]).reshape(
        int(math.sqrt(len(r_values))),
        int(math.sqrt(len(t_values))))
    plt.contour(r_values, t_values, npz)
    plt.show()


if __name__ == "__main__":
    R_min = 1
    R_max = 5
    N = 50

    r_step = (R_max - R_min)/(N+1)
    t_step = (2 * math.pi - 0)/N
    r_values = seq(np.arange(R_min+r_step, R_max, r_step).tolist()) \
        .map(lambda r : [r]*N) \
        .reduce(lambda a,b: a+b) \
        .to_list()
    print(len(r_values))
    t_values = list(reduce(
        (lambda a,b:a+b),
        [np.arange(0, 2*math.pi, t_step).tolist()]*N))
    print(len(t_values))

    z_values = list(map(psi_exact, r_values, t_values))

    polar_plot(r_values, t_values, z_values)
    contour_plot(r_values, t_values, z_values)
