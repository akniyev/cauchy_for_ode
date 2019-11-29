from laguerre_functions import *


class CauchySolver:
    def __init__(self, n=10, alpha=0, epsilon=1e-12):
        super()
        self.alpha = alpha
        self.n = n
        self.epsilon = epsilon

    def roots(self):
        return find_all_roots_of_laguerre(self.n, self.alpha, self.epsilon)

    def laguerre_on_grid(self, a, b, density=2000):
        xs = [a + (b - a)/density * i for i in range(density+1)]
        ys = [laguerre(self.n, self.alpha, x) for x in xs]
        return xs, ys
