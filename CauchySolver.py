from typing import List

from laguerre_functions import *


class CauchySolver:
    def __init__(self, n=10, alpha=0, n_part=15, a=1, b=1, epsilon=1e-12):
        super()
        self.alpha = alpha
        self.n = n
        self.n_part = n_part
        self.a = a
        self.b = b
        self.epsilon = epsilon
        self.iteration_number = 0
        self.current_cs = []

        self.y0 = 1

    def set_parameters(self, n=None, alpha=None, n_part=None, a=None, b=None, epsilon=None):
        self.n = n if n is not None else self.n
        self.n_part = n_part if n_part is not None else self.n_part
        self.alpha = alpha if alpha is not None else self.alpha
        self.a = a if a is not None else self.a
        self.b = b if b is not None else self.b
        self.epsilon = epsilon if epsilon is not None else self.epsilon

    def roots(self) -> List[float]:
        return find_all_roots_of_laguerre(self.n, self.alpha, self.epsilon)

    def laguerre_on_grid(self, a, b, density=2000):
        xs = [a + (b - a)/density * i for i in range(density+1)]
        ys = [laguerre(self.n, self.alpha, x) for x in xs]
        return xs, ys

    def reset_iterations(self):
        self.current_cs = [1/(i+1) for i in range(self.n_part+1)]
        self.iteration_number = 0

    def next_iterations(self):
        new_cs = perform_iteration_on_cs(self.current_cs, self.a, self.b, self.alpha, self.n_part, self.n)
        self.iteration_number += 1
        self.current_cs = new_cs

    def get_solution(self, density):
        a = self.a
        b = self.b
        n = self.n
        n_part = self.n_part

        xs = [1/density for i in range(density)]
        ys = [find_solution(-math.log(1-x)/a, n, a, b, self.alpha, self.current_cs, n_part, self.y0) for x in xs]

        return xs, ys