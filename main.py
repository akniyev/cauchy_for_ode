laguerre_cache = dict()
epsilon = 1e-10
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
    mid_point = (a + b) / 2
    fa = laguerre(k, alpha, a)
    fb = laguerre(k, alpha, b)
    fc = laguerre(k, alpha, mid_point)

    if abs(fa - fb) < epsilon:
        return (fa + fb) / 2

    if fa * fb > 0:
        raise Exception("The function should have different signs at the ends")

    if abs(fc) < epsilon:
        return mid_point
    elif fa * fc < 0:
        return find_root_of_laguerre(k, alpha, a, mid_point)
    else:
        return find_root_of_laguerre(k, alpha, mid_point, b)


n = int(input("Enter the order: "))

max_bound = 4*n + 2

segments = {1: [(0, max_bound)]}
roots = dict()

for k in range(1, n+1):
    layer_roots = []
    for s in segments[k]:
        left_end = s[0]
        right_end = s[1]
        segment_root = find_root_of_laguerre(k, alpha, left_end, right_end)
        layer_roots.append(segment_root)

    next_segments = []
    next_segments.append((0, layer_roots[0]))
    for i in range(len(layer_roots) - 1):
        next_segments.append((layer_roots[i], layer_roots[i+1]))
    next_segments.append((layer_roots[-1], max_bound))

    segments[k+1] = next_segments
    roots[k] = layer_roots

print(roots[n])
