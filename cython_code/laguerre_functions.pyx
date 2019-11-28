# cython: language_level=3, boundscheck=False


cdef double laguerre_iterative_cython(int k, double alpha, double x):
    if k <= 0:
        return 1

    cdef double minus_2 = 1
    cdef double minus_1 = -x + alpha + 1
    cdef int i
    cdef double a
    cdef double b
    cdef double current

    for i in range(2, k+1):
        a = (2.0 * i - 1 + alpha - x) / i
        b = (i + alpha - 1.0) / i
        current = a * minus_1 - b * minus_2
        minus_2 = minus_1
        minus_1 = current

    return minus_1


cdef double find_root_of_laguerre_cython(int k, double alpha, double a, double b, double epsilon):
    """
    Finds the root of the Laguerre polynomials of order k on the segment [a,b].
    The root must be single.
    """
    cdef double mid_point = (a + b) / 2
    cdef double fa = laguerre_iterative_cython(k, alpha, a)
    cdef double fc = laguerre_iterative_cython(k, alpha, mid_point)

    if abs(a - b) < epsilon or abs(fc) < epsilon:
        return mid_point

    if fa * fc < 0:
        return find_root_of_laguerre_cython(k, alpha, a, mid_point, epsilon)
    else:
        return find_root_of_laguerre_cython(k, alpha, mid_point, b, epsilon)


def laguerre(k, alpha, x):
    return laguerre_iterative_cython(k, alpha, x)


def find_root_of_laguerre(k, alpha, a, b, epsilon):
    return find_root_of_laguerre_cython(k, alpha, a, b, epsilon)
