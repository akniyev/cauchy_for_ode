import math
import matplotlib

epsilon = 1e-12
alpha = 0


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

    a = (2 * k - 1 + alpha - x) / k
    b = (k + alpha - 1) / k

    return a * laguerre(k - 1, alpha, x) - b * laguerre(k - 2, alpha, x)


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
    :return: the array containing all the roots
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
def laguerre_derivative(k, alpha, x):
    return -laguerre(k - 1, alpha + 1, x)


@cached
def root_to_weight(root, n):
    return 1.0 / root / ((laguerre(n - 1, alpha + 1, root)) ** 2)


@cached
def sobolev_laguerre(tau, b, j):
    return math.sqrt(b) * tau * laguerre(j - 1, 1, b * tau) / j


def etta(x):
    return 1


def f(x, y):
    return 1


@cached
def g(tau, k, a, b):
    first_arg = 1 - math.exp(-a * tau)
    second_arg_sum_c_sobolev_laguerre = etta(0) + sum(
        [c[j] * sobolev_laguerre(tau, b, j + 1) for j in range(n_part + 1)])
    return f(first_arg, second_arg_sum_c_sobolev_laguerre) * laguerre(k, 0, tau * b) * math.exp((1 - a - b) * tau)


# n = int(input("Enter the order: "))
n = 10
n_part = 10

c = [1 / n_part for i in range(n_part)]

a = 1
b = 1
k = 2

roots = find_all_roots_of_laguerre(n, alpha)
weights = [root_to_weight(root, n) for root in roots]

result = sum([(roots[i] ** 2) * weights[i] for i in range(n)])
print(result)

# result = a * sum([weights[i] * g(roots[i], k, a, b) for i in range(n)])


print(roots)
print(weights)



