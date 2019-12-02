import cython_code.laguerre_functions as lf
import math


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


def laguerre(k, alpha, x):
    return lf.laguerre(k, alpha, x)


def find_root_of_laguerre(k, alpha, a, b, epsilon):
    """
    Finds the root of the Laguerre polynomials of order k on the segment [a,b].
    The root must be single.
    """
    return lf.find_root_of_laguerre(k, alpha, a, b, epsilon)


def find_all_roots_of_laguerre(k, alpha, epsilon = 1e-12):
    """
    Finds all k roots of the Laguerre polynomials of order k
    :param k: the order of the polynomials
    :param alpha: a parameter of the Laguerre polynomials
    :return: the array containing all the k roots roots[0], ...,roots[k-1]
    """
    max_bound = 4 * k + 2

    segments = list(range(k+2))
    segments[1] = [(0, max_bound)]
    roots = list(range(k+2))

    for j in range(1, k + 1):
        layer_roots = []
        for s in segments[j]:
            left_end = s[0]
            right_end = s[1]
            segment_root = find_root_of_laguerre(j, alpha, left_end, right_end, epsilon)
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
    return 1.0 / root / ((laguerre(k - 1, 1, root)) ** 2)


@cached
def all_weights(_k, _alpha):
    _roots = find_all_roots_of_laguerre(_k, _alpha)
    return [root_to_weight(root, _k, _alpha) for root in _roots]


@cached
def sobolev_laguerre(tau, b, j):
    return math.sqrt(b) * tau * laguerre(j - 1, 1, b * tau) / j


def g(tau, _k, _a, _b, _c, _alpha, n_part, f, y0):
    first_arg = 1 - math.exp(-_a * tau)
    s0 = sum([_c[j] * sobolev_laguerre(tau, _b, j + 1) for j in range(n_part + 1)])
    second_arg_sum_c_sobolev_laguerre = y0 + s0
    s1 = math.sqrt(_b) * laguerre(_k, 0, tau * _b) * math.exp((1 - _a - _b) * tau)
    return f(first_arg, second_arg_sum_c_sobolev_laguerre) * s1


def perform_iteration_on_c(c0, _a, _b, _k, _alpha, n_part, _n):
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


def find_solution(t, _k, _a, _b, _alpha, c, n_part, _n, y0):
    return y0 + sum([c[k] * sobolev_laguerre(t, _b, k+1) for k in range(n_part+1)])