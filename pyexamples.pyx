cdef extern from "examples.h":
    void hello(const char *name)

cdef extern from "examples.h":
    double d_abs(double x)

def py_hello(name: bytes) -> None:
    hello(name)

def py_d_abs(x: double) -> double:
    return d_abs(x)