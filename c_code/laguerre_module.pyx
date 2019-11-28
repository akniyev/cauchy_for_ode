cdef extern from "laguerre.h":
    double laguerre_iterative(int k, double alpha, double x)


def py_laguerre_iterative(k: int, alpha: double, x: double) -> double:
    return laguerre_iterative(k, alpha, x)
