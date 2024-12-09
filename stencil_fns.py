import sys
import math
import numpy as np
from functional import seq
from functools import reduce

def center(r,dr,dt):
    return -2 * (1/pow(dr,2) + (1/(pow(r,2) * pow(dt,2))))

def upper(r,dr):
    return (1/pow(dr,2))+(1/(2*r*dr))

def lower(r,dr):
    return (1/pow(dr,2))-(1/(2*r*dr))

def left_right(r,dt):
    return (1/(pow(r,2) * pow(dt,2)))

def stencil_check(dr,dt,r):
    print("center @r="+str(r)+": "+str(center(r,dr,dt)))
    print("upper  @r="+str(r)+": "+str(upper(r,dr)))
    print("lower  @r="+str(r)+": "+str(lower(r,dr)))
    print("left   @r="+str(r)+": "+str(left_right(r,dt)))
    print("right  @r="+str(r)+": "+str(left_right(r,dt)))

def outer_bc(R_max, t):
    return R_max * math.sin(t)

def d_range(start, end, N, st_ex=False):
    step = (end-start)/(N+1 if st_ex else N)
    np_range = np.arange(
        (start+step if st_ex else start),
        end,
        step
    )
    return np_range[0:N].tolist()

def r_values(R_min, R_max, N):
    return seq(d_range(R_min, R_max, N, st_ex=True)) \
    .map(lambda r : [r]*N) \
    .reduce(lambda a,b: a+b) \
    .to_list()

def t_values(N):
    return list(reduce(
        (lambda a,b:a+b),
        [d_range(0, 2*math.pi, N)]*N))


if __name__ == "__main__":
    dr = 0.8
    dt=1.5708
    R_max = 5
    match sys.argv[1]:
        case "stencil":
            r = float(sys.argv[2])
            print("r = " + str(r))
            stencil_check(dr,dt,r)
        case "bc":
            t = float(sys.argv[2])
            print("t = " + str(t))
            print("outer_bc @t=" + str(t)+": "+str(outer_bc(R_max, t)))
