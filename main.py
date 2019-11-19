import math
import matplotlib.pyplot as plt
import numpy as np

epsilon = 1e-12
# alpha = 0

# START: Initial conditions
y0 = 1


def f(x, y):
    return x * y


# END: Initial conditions

def cached(func):
    """
    A decorator which caches values for pure functions
    """
    cache = {}

    def wrapper(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return wrapper


@cached
def laguerre(k, alpha, x):
    if k == 0:
        return 1
    elif k == 1:
        return -x + alpha + 1

    _a = (2 * k - 1 + alpha - x) / k
    _b = (k + alpha - 1) / k

    return _a * laguerre(k - 1, alpha, x) - _b * laguerre(k - 2, alpha, x)


@cached
def find_root_of_laguerre(k, alpha, a, b):
    """
    Finds the root of the Laguerre polynomials of order k on the segment [a,b].
    The root must be single.
    """
    mid_point = (a + b) / 2
    fa = laguerre(k, alpha, a)
    fb = laguerre(k, alpha, b)
    fc = laguerre(k, alpha, mid_point)

    if abs(a - b) < epsilon:
        return (a + b) / 2

    if fa * fb > 0:
        raise Exception("The function should have different signs at the ends")

    if abs(fc) < epsilon:
        return mid_point

    elif fa * fc < 0:
        return find_root_of_laguerre(k, alpha, a, mid_point)
    else:
        return find_root_of_laguerre(k, alpha, mid_point, b)


@cached
def find_all_roots_of_laguerre(k, alpha):
    """
    Finds all k roots of the Laguerre polynomials of order k
    :param k: the order of the polynomials
    :param alpha: a parameter of the Laguerre polynomials
    :return: the array containing all the k roots roots[0], ...,roots[k-1]
    """
    max_bound = 4 * k + 2

    segments = {1: [(0, max_bound)]}
    roots = dict()

    for j in range(1, k + 1):
        layer_roots = []
        for s in segments[j]:
            left_end = s[0]
            right_end = s[1]
            segment_root = find_root_of_laguerre(j, alpha, left_end, right_end)
            layer_roots.append(segment_root)

        next_segments = [(0, layer_roots[0])]
        for i in range(len(layer_roots) - 1):
            next_segments.append((layer_roots[i], layer_roots[i + 1]))
        next_segments.append((layer_roots[-1], max_bound))

        segments[j + 1] = next_segments
        roots[j] = layer_roots

    return roots[k]


@cached
def root_to_weight(root, k, alpha):
    """
    Given one root of the Laguerre polynomials of order k, calculates the corresponding weight
    :param root:
    :param k: int
    :return:
    """
    return 1.0 / root / ((laguerre(k - 1, alpha + 1, root)) ** 2)


@cached
def sobolev_laguerre(tau, b, j):
    return math.sqrt(b) * tau * laguerre(j - 1, 1, b * tau) / j


@cached
def all_weights(_k, _alpha):
    _roots = find_all_roots_of_laguerre(_k, _alpha)
    return [root_to_weight(root, _k, _alpha) for root in _roots]


def g(tau, _k, _a, _b, _c, _alpha, n_part):
    first_arg = 1 - math.exp(-_a * tau)
    second_arg_sum_c_sobolev_laguerre = y0 + sum([_c[j] * sobolev_laguerre(tau, _b, j + 1) for j in range(n_part + 1)])
    return f(first_arg, second_arg_sum_c_sobolev_laguerre) * laguerre(_k, 0, tau * _b) * math.exp((1 - _a - _b) * tau)


def perform_iteration_on_c(c0, _a, _b, _k, _alpha, n_part, _n):
    """
    :param c0: a _k+1 array
    :param _a:
    :param _b:
    :param _k:
    :param _alpha:
    :param n_part:
    :return:
    """
    weights = all_weights(_n, _alpha)
    _roots = find_all_roots_of_laguerre(_n, _alpha)
    gs = [g(root, _k, _a, _b, c0, _alpha, n_part) for root in _roots]

    return _a * sum([weights[i] * gs[i] for i in range(_n)])


def perform_iteration_on_cs(c0, _a, _b, _alpha, n_part, _n):
    return [perform_iteration_on_c(c0, _a, _b, k, _alpha, n_part, _n) for k in range(n_part+1)]


def distance(c0, c1):
    _n = len(c0)
    return math.sqrt(sum([(c0[i] - c1[i]) ** 2 for i in range(_n)]))


def find_cs(_k, _a, _b, _alpha, threshold, n_part, _n):
    c0 = [1/(i+1) for i in range(n_part+1)]
    c1 = perform_iteration_on_cs(c0, _a, _b, _alpha, n_part, _n)

    while distance(c0, c1) > threshold:
        c0 = c1
        c1 = perform_iteration_on_cs(c0, _a, _b, _alpha, n_part, _n)

    return c1


def find_solution(t, _k, _a, _b, _alpha, threshold, n_part, _n):
    c = find_cs(_k, _a, _b, _alpha, threshold, n_part, _n)
    return y0 + sum([c[k] * sobolev_laguerre(t, _b, k+1) for k in range(n_part+1)])


n_dense = 2000
c_threshold = 1e-3
n = 20
alpha = 0
n_part = 15
a = 1
b = 1

ts = [i/n_dense for i in range(n_dense)]
ys = [find_solution(-math.log(1-t)/a, n, a, b, alpha, c_threshold, n_part, n) for t in ts]
solution_ys = [math.exp((t ** 2)/2) for t in ts]


plt.plot(ts, ys)
plt.plot(ts, solution_ys)
plt.show()
