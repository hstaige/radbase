import inspect
import time
from abc import ABC, abstractmethod
from functools import lru_cache
from typing import Callable, Literal, TypeAlias

import jax
import jax.numpy as jnp
import numpy as np
import numpy.typing as npt
from jax import Array
from jax.typing import ArrayLike
from pint import UnitRegistry
from scipy.integrate import quad
from sympy import Symbol, lambdify, parse_expr

jax.config.update('jax_enable_x64', True)


def get_moment_unit(moment: str, ureg: UnitRegistry):
    moment_expr = compile_moment_expr(moment) if isinstance(moment, str) else moment
    if 'theta' in inspect.signature(moment_expr).parameters:
        return moment_expr(r=ureg('fm'), theta=ureg(''))
    else:
        return moment_expr(r=ureg('fm'))


@lru_cache(maxsize=1000)
def compile_moment_expr(expression: str) -> Callable[[ArrayLike], Array] | Callable[[ArrayLike, ArrayLike], Array]:
    r, theta = Symbol('r', nonnegative=True), Symbol(r'theta')
    expression = parse_expr(expression, {'r': r, 'theta': theta})
    symbols = expression.atoms(Symbol)

    if unknown := (symbols - {r, theta}):
        raise ValueError(f'Invalid symbols in {expression}: found {unknown}, symbols must only include {r} or {theta}.')

    symbols = sorted(symbols, key=lambda s: 0 if s == r else 1)  # preserve r, theta ordering
    return lambdify(symbols, expression, 'jax')


def trapezoid_weights(num_points: int, deltax, datatype):
    w = jnp.ones((num_points,), dtype=datatype)
    w = w.at[0].set(0.5)
    w = w.at[-1].set(0.5)
    return w * deltax


MomentExpr: TypeAlias = str | Callable[[ArrayLike], Array]


class ChargeDistribution(ABC):
    name: str
    parameter_names: tuple[str, ...]
    test_r: npt.NDArray

    _fcn_type: Literal['spherical', 'deformed']
    _use_radial_segments: bool = False

    def __init__(self, name: str, parameter_names: tuple[str, ...],
                 test_r: npt.NDArray):
        self.name = name
        self.parameter_names = parameter_names
        self.test_r = test_r
        # moment str / function id, r_cut, n_r_points, n_theta_points
        self._quadrature_grid_cache = {}
        self._moment_vg_cache: dict[tuple[str | int, float, int, int], Callable[[dict[str, float]], Array]] = {}

    def get_radial_segments(self, theta, pars, r_cut: float, dtype=jnp.float64):
        """
        Returns bounds of regions where charge distribution is smooth. Integration is done independently in each region.

        Should be overridden by charge distributions with discontinuities.

        Example: Uniform charge distribution with varying cutoff radius

        lower = jnp.array([0.0], dtype=dtype)
        upper = jnp.array([pars['R'], dtype=dtype)
        return lower, upper
        """
        lower = jnp.array([0.0], dtype=dtype)
        upper = jnp.array([r_cut], dtype=dtype)
        return lower, upper

    def density_on_segments(self, r, theta, pars: dict[str, float], segment_index: int):
        return self(r, theta, pars)

    @abstractmethod
    def __call__(self, r, theta, pars) -> Array:
        pass

    def _get_moment_value_and_gradient_fcn(self,
                                           moment: MomentExpr,
                                           r_cut: float = 50.0,
                                           n_r_points: int = 1024,
                                           n_theta_points: int = 128,
                                           dtype=jnp.float64):

        moment_key = moment if isinstance(moment, str) else id(moment)
        key = (moment_key, r_cut, n_r_points, n_theta_points)
        if key in self._moment_vg_cache:
            return self._moment_vg_cache[key]

        moment_expr = compile_moment_expr(moment) if isinstance(moment, str) else moment

        sig = inspect.signature(moment_expr)
        if unknown := set(sig.parameters) - {'r', 'theta'}:
            raise ValueError(f'Moment function has unknown parameters {unknown}, must only contain r and theta.')
        elif 'r' not in sig.parameters:
            raise ValueError(f'Moment function must contain r in the parameters. Only contains {set(sig.parameters)}.')

        # moment_uses_theta = True if 'theta' in sig.parameters else False

        if 'theta' not in sig.parameters and self._fcn_type == 'spherical':
            # noinspection PyTypeChecker
            r: Array = jnp.linspace(0.0, r_cut, n_r_points, dtype=dtype)

            # Assuming linearly spaced trapezoidal pieces, could multiply by dx but washes out in norm
            weights = trapezoid_weights(n_r_points, r_cut / (n_r_points - 1), dtype)
            weights = weights * r * r
            moment_weights = weights * moment_expr(r)

            def core(pars):
                grid = self(r, 0, pars=pars)
                norm = jnp.sum(grid * weights)
                return jnp.sum(grid * moment_weights) / norm

            res = jax.jit(jax.value_and_grad(core))
        else:
            # noinspection PyTypeChecker
            r: Array = jnp.linspace(0.0, r_cut, n_r_points, dtype=dtype)
            theta = jnp.linspace(0, jnp.pi, n_theta_points)
            r, theta = jnp.meshgrid(r, theta, indexing='ij')

            r_weights = trapezoid_weights(n_r_points, r_cut / (n_r_points - 1), dtype).reshape((1, n_r_points))
            theta_weights = trapezoid_weights(n_theta_points, jnp.pi / (n_theta_points - 1), dtype).reshape(
                (1, n_theta_points))
            weights = r_weights.T @ theta_weights
            weights = weights * r * r * jnp.sin(theta)

            moment_weights = weights * moment_expr(r, theta) if 'theta' in sig.parameters else weights * moment_expr(r)

            def core(pars):
                grid = self(r, theta, pars=pars)  # spherical functions allow theta but ignore it
                norm = jnp.sum(grid * weights)
                return jnp.sum(grid * moment_weights) / norm

            res = jax.jit(jax.value_and_grad(core))

        self._moment_vg_cache[key] = res
        return res

    def calculate_moment_gradient(self, moment, pars: dict[str, float], r_cut: float = 50.0, n_r_points: int = 1024,
                                  n_theta_points: int = 128):
        return self._get_moment_value_and_gradient_fcn(moment, r_cut=r_cut, n_r_points=n_r_points,
                                                       n_theta_points=n_theta_points)(pars)

    def calculate_moment_gradient_scipy(self, moment, pars: dict[str, float], r_cut: float = 50.0):

        moment_expr = compile_moment_expr(moment) if isinstance(moment, str) else moment

        def calculate_gradient_scipy(parameters, weighting: Literal['r**2', ''] = 'r**2'):
            value, error, *_ = quad(lambda r: moment_expr(r) * r ** 2 * self.__call__(r, None, pars=parameters),
                                    0, r_cut, epsrel=1e-5)
            norm, norm_error, *_ = quad(lambda r: r ** 2 * self.__call__(r, None, pars=parameters),
                                        0, r_cut, epsrel=1e-5)
            return value

        moment_value = calculate_gradient_scipy(pars) / calculate_gradient_scipy(pars, weighting='')
        return moment_value


class TwoParameterFermi(ChargeDistribution):

    def __init__(self, name: str = '2pf', parameter_names: tuple[str, ...] = ('c', 't'),
                 test_r: npt.NDArray | None = None):
        self._fcn_type = 'spherical'
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def __call__(self, r: Array, _, pars: dict[str, float]) -> Array:
        c, t = pars['c'], pars['t']
        return 4 * jnp.pi * 1 / (1 + jnp.exp(4 * jnp.log(3) * (r - c) / t))


class ThreeParameterFermi(ChargeDistribution):

    def __init__(self, name: str = '3pf', parameter_names: tuple[str, ...] = ('c', 't', 'beta2'),
                 test_r: npt.NDArray | None = None):
        self._fcn_type = 'deformed'
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def __call__(self, r: Array, theta: Array, pars: dict[str, float]) -> Array:
        c, t, beta2 = pars['c'], pars['t'], pars['beta2']
        return 4 * jnp.pi * 1 / (1 + jnp.exp(4 * jnp.log(3) * (r - c * (1 + beta2 * self.sph_harm_20(theta))) / t))

    def sph_harm_20(self, theta):
        return 1 / 4 * jnp.sqrt(5 / jnp.pi) * (3 * jnp.cos(theta) ** 2 - 1)


class ParabolicFermi(ChargeDistribution):

    def __init__(self, name: str = 'parabolic_fermi', parameter_names: tuple[str, str, str] = ('c', 't', 'w'),
                 test_r: npt.NDArray | None = None):
        self._fcn_type = 'spherical'
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def __call__(self, r: Array, _, pars: dict[str, float]) -> Array:
        c, t, w = pars['c'], pars['t'], pars['w']
        return 4 * jnp.pi * (1 + w * r ** 2 / c ** 2) / (1 + jnp.exp(4 * jnp.log(3) * (r - c) / t))


class Uniform(ChargeDistribution):

    def __init__(self, name: str = 'uniform', parameter_names: tuple[str] = ('R',),
                 test_r: npt.NDArray | None = None):
        self._fcn_type = 'spherical'
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def __call__(self, r: Array, _, pars: dict[str, float]) -> Array:
        r_cutoff = pars['R']
        return -jnp.tanh((r - r_cutoff) * 100) + 1


cds = {cd.name: cd for cd in [TwoParameterFermi(), ThreeParameterFermi(), ParabolicFermi()]}

if __name__ == '__main__':
    num_runs = 10000

    twopf = TwoParameterFermi()
    test_par = {'c': 5.0, 't': 2.3}

    uniform = Uniform()
    print(uniform.calculate_moment_gradient('r**2', pars={'R': 5.0}))

    uniform = Uniform()
    print(uniform.calculate_moment_gradient('r**2', pars={'R': 5.0}, n_r_points=2048))

    print(twopf.calculate_moment_gradient_scipy('r**2', test_par))
    start = time.time()
    for _ in range(500):
        twopf.calculate_moment_gradient_scipy('r**2', test_par)
    post_twopf_scipy = time.time()
    print(f'Took {(time_taken := (post_twopf_scipy - start)): 0.3g} s for {num_runs} integrals,'
          f' or {time_taken / num_runs: 0.3g} s per integral.')

    start = time.time()
    vg_func = jax.block_until_ready(twopf._get_moment_value_and_gradient_fcn('r**2'))
    print(f'{time.time() - start: 0.3g} s to compile the moment function.')

    start = time.time()
    for _ in range(num_runs):
        value, gradient = vg_func(test_par)
    jax.block_until_ready(gradient)
    post_twopf = time.time()
    print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
          f' or {time_taken / num_runs: 0.3g} s per integral.')
    print(value, gradient)

    print('Same but more trapezoids.')
    start = time.time()
    vg_func = jax.block_until_ready(twopf._get_moment_value_and_gradient_fcn('r**2', n_r_points=4096))
    print(f'{time.time() - start: 0.3g} s to compile the moment function.')

    num_runs = 10000
    start = time.time()
    for _ in range(num_runs):
        value, gradient = vg_func(test_par)
    jax.block_until_ready(gradient)
    post_twopf = time.time()
    print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
          f' or {time_taken / num_runs: 0.3g} s per integral.')
    print(value, gradient)

    print('\nAnd now for the deformed case:')
    threepf = ThreeParameterFermi()
    test_par = {'c': 5.0, 't': 2.3, 'beta2': 0.3}

    start = time.time()
    vg_func = jax.block_until_ready(threepf._get_moment_value_and_gradient_fcn('r**2'))
    print(f'{time.time() - start: 0.3g} s to compile the moment function.')

    start = time.time()
    for _ in range(num_runs):
        value, gradient = vg_func(test_par)
    jax.block_until_ready(gradient)
    post_twopf = time.time()
    print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
          f' or {time_taken / num_runs: 0.3g} s per integral.')
    print(value, gradient)
