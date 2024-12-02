import ctypes
import pathlib
from numpy.ctypeslib import ndpointer

if __name__ == "__main__":
    libname=pathlib.Path().absolute() / "laplace-operator-2"
    c_lib = ctypes.CDLL(libname)

    polar_laplace = c_lib.polar_laplace
    #polar_laplace.restype = ctypes.POINTER(ctypes.c_double * 8)
    polar_laplace.restype = ndpointer(dtype=ctypes.c_double, shape=(8,))
    #polar_laplace.restype = ctypes.c_void_p

    answer = polar_laplace(1, 5, 3)
#    print("here")
#    d_arr = ctypes.cast(answer, ctypes.POINTER(ctypes.c_double))
#    c_lib.freeMem(answer)
