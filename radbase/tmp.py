from functools import partial
from timeit import repeat
from typing import Callable

import numpy as np
import numpy.typing as npt
from scipy.integrate import trapezoid
from sympy import Symbol, lambdify, parse_expr


def f2p(r, c, t):
    return 1 / (1 + np.exp(-(r - c) * t))


def compile_expr(expression):
    expression = parse_expr(expression)
    symbols = tuple(expression.atoms(Symbol))
    return lambdify(symbols, expression, 'jax')


def integrate_trap_sympy(moment_function, c, t, r_cut, n_points=100):

    x = np.linspace(0, r_cut, n_points)
    y = f2p(x, c, t) * moment_function(r=x)
    return trapezoid(y, x)


def integrate_trap_numpy(r_pow, alpha, c, t, r_cut, n_points=100):
    x = np.linspace(0, r_cut, n_points)
    y = f2p(x, c, t) * x ** r_pow * np.exp(-alpha * x)
    return trapezoid(y, x)


def eval_sympy(r: npt.NDArray, compiled: Callable[[npt.NDArray], npt.NDArray]) -> npt.NDArray:
    return compiled(r)


def get_time(func, *args, number=100, repeat_num=5):
    rep_data = np.array(repeat(partial(func, *args), number=number, repeat=repeat_num))
    return rep_data.min(), rep_data.mean()


if __name__ == '__main__':
    moment_expr = compile_expr('r**2')
    print(get_time(eval_sympy, np.linspace(0, 100, 100), moment_expr))

    print(get_time(integrate_trap_sympy, compile_expr('r**2'), 5.2, 2.3, 30, 10000))
    print(get_time(integrate_trap_numpy, 2, 0, 5.2, 2.3, 30, 10000))
