import time

laguerre_cache = dict()
epsilon = 1e-12
alpha = 0


def laguerre(k, alpha, x):
    if (k, alpha, x) in laguerre_cache:
        return laguerre_cache[(k, alpha, x)]
    if k == 0:
        return 1
    elif k == 1:
        return -x + alpha + 1

    a = (2*k - 1 + alpha - x) / k
    b = (k + alpha - 1) / k

    result = a * laguerre(k - 1, alpha, x) - b * laguerre(k - 2, alpha, x)
    laguerre_cache[(k, alpha, x)] = result
    return result


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


def laguerre_derivative(k, alpha, x):
    return -laguerre(k-1, alpha+1, x)


n = int(input("Enter the order: "))

roots = find_all_roots_of_laguerre(n, alpha)






