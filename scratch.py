import math
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from functional import seq #TODO: Delete
from functools import reduce #TODO: Delete

def d_range(start, end, N, st_ex=False):
    step = (end-start)/(N+1 if st_ex else N)
    np_range = np.arange(
        (start+step if st_ex else start),
        end,
        step
    )
    return np_range[0:N].tolist()

def sample_data():
    R_min = 1
    R_max = 5
    N = 5
    r_values = seq(d_range(R_min, R_max, N, st_ex=True)) \
        .map(lambda r : [r]*N) \
        .reduce(lambda a,b: a+b) \
        .to_list()
    t_values = list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))
    z_values = list(map(
        (lambda r,t : r * math.sin(t)),
        r_values, t_values
    ))
    return {"r" : r_values, "t" : t_values, "z" : z_values}


def contour_plot_old(r_values, t_values, z_values):
    print("r_values: " + str(r_values))
    print("t_values: " + str(t_values))
    print("\n\n")
    unique_r = list(set(r_values))
    unique_t = list(set(t_values))
    #unique_r += [unique_r[0]]
    #unique_t += [unique_t[0]]
    print("unique_r: " + str(unique_r))
    print("unique_t: " + str(unique_t))

    r,t = np.meshgrid(unique_r, unique_t)
    z = np.array(z_values).reshape(len(unique_r), len(unique_t)) \
                          .transpose()
    print(z)
    #z = np.append(z, [z[0,:]], axis=0)
    #print(z)
    #first_col = list(map((lambda i : [i]), z[:,0]))
    #print(first_col)
    #print("[[0]]*len(unique_r): " + str([[0]]*len(unique_r)))
    #z = np.append(z, [[0]]*len(unique_r), axis=1)
    #z = np.append(z, first_col, axis=1)
    #print(z)

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_ylim(0,5)
    plt.plot([0,0], [0,2], 'ro')
    ax.contour(t, r, z)
    circle = pl.Circle((0.0, 0.0), 1.0, transform=ax.transData._b, color="red", alpha=0.4)
    ax.add_artist(circle)
    plt.show()

def contour_plot(r_values, t_values, z_values):
    unique_r = list(set(r_values))
    unique_t = list(set(t_values))
    p2c = (lambda r,t : (r*math.cos(t), r*math.sin(t)))
    x_values, y_values = zip(*list(map(p2c, unique_r, unique_t)))
    print(len(x_values))
    print(len(y_values))
    x,y = np.meshgrid(x_values, y_values)

    z = np.array(z_values).reshape(len(unique_r), len(unique_t)).transpose()
    print(z)

    fig, ax = plt.subplots()
    ax.set_xlim(-5,5)
    ax.set_ylim(-5,5)
    ax.contour(x,y,z)
    plt.show()

