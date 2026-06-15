import inspect
import re
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
from sympy import Symbol, lambdify, parse_expr

jax.config.update('jax_enable_x64', True)


def get_moment_unit(moment: str | Callable, ureg: UnitRegistry):
    moment_expr = compile_moment_expr(moment) if isinstance(moment, str) else moment
    if 'theta' in inspect.signature(moment_expr).parameters:
        return moment_expr(r=ureg('fm'), theta=ureg(''))
    else:
        return moment_expr(r=ureg('fm'))


@jax.jit
def r2(r):
    return r * r


@jax.jit
def rk(r, k):
    return r ** k


@jax.jit
def rka(r, k, alpha):
    return r ** k * jnp.exp(-alpha * r)


# positive decimal number is \d*[.,]?\d*
rk_match = re.compile(r'^r\*\*(?P<k>-?d*\.?\d*)$')
rka_match = re.compile(r'^r\*\*(?P<k>-?d*\.?\d*)\*exp(-(?P<alpha>\d*\.?\d*)\*r)$')


@lru_cache(maxsize=1000)
def compile_moment_expr(expression: str) -> Callable[[ArrayLike], Array] | Callable[[ArrayLike, ArrayLike], Array]:
    if match := rk_match.match(expression):
        k = float(match.group('k'))
        if k == 2.0:
            return r2
        return jax.tree_util.Partial(rk, k=k)
    elif match := rka_match.match(expression):
        k, alpha = float(match.group('k')), float(match.group('alpha'))
        return jax.tree_util.Partial(rka, k=k, alpha=alpha)

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


def _reshape_segment_bound(x, angular: bool, dtype):
    """
    Normalize lower/upper segment-bound arrays to broadcastable shapes.

    If non-angular, target shape is (n_segments, 1)

    If angular but angle independent bounds, target_shape is (n_segment, 1, 1)

    If angular but angle dependent bounds, target_shape is (n_segment, 1, n_theta)
    """
    x = jnp.asarray(x, dtype=dtype)

    if angular:
        if x.ndim == 0:
            return x.reshape((1, 1, 1))

        if x.ndim == 1:
            # theta-independent bounds: shape (n_segments,)
            return x.reshape((-1, 1, 1))

        if x.ndim == 2:
            # theta-dependent bounds: shape (n_segments, n_theta)
            return x[:, None, :]

        raise ValueError(
            "Angular segment bounds must be scalar, shape (n_segments,), "
            "or shape (n_segments, n_theta)."
        )

    else:
        if x.ndim == 0:
            return x.reshape((1, 1))

        if x.ndim == 1:
            return x[:, None]

        if x.ndim == 2 and x.shape[1] == 1:
            return x

        raise ValueError(
            "Spherical segment bounds must be scalar, shape (n_segments,), "
            "or shape (n_segments, 1)."
        )


MomentExpr: TypeAlias = str | Callable[[ArrayLike], Array]


class ChargeDistribution(ABC):
    name: str
    parameter_names: tuple[str, ...]
    test_r: npt.NDArray

    _fcn_type: Literal['spherical', 'deformed']
    _use_radial_segments: bool = False
    _n_radial_segments: int = 1

    def __init__(self, name: str, parameter_names: tuple[str, ...],
                 test_r: npt.NDArray):
        self.name = name
        self.parameter_names = parameter_names
        self.test_r = test_r
        # moment str / function id, r_cut, n_r_points, n_theta_points
        self._quadrature_grid_cache = {}
        self._moment_vg_cache: dict[tuple[str | int, float, int, int, str], Callable[[dict[str, float]], Array]] = {}
        self._fixed_vg_cache: dict[tuple[bool, bool, float, int, int, str], tuple[Array, Array, Callable]] = {}

    def radial_segments(self, theta, pars, r_cut: float, dtype):
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

    def density_on_segments(self, r, theta, pars: dict[str, float], segment_index):
        return self(r, theta, pars)

    @abstractmethod
    def __call__(self, r, theta, pars) -> Array:
        pass

    def _get_segment_reference_grid(self, angular: bool, n_r_points: int = 256, n_theta_points: int = 64,
                                    dtype=jnp.float64):
        """
        If angular is False, grid will have shape (num_segments, n_r_points).
        If angular is True, grid will have shape (num_segments, n_r_points, n_theta_points)

        Provides the base, unscaled segments that will later be converted to the actual r and theta regions.
        """

        key = (angular, n_r_points, n_theta_points, str(dtype))

        if key in self._quadrature_grid_cache:
            return self._quadrature_grid_cache[key]

        n_segments = self._n_radial_segments

        if not angular:
            u_np, u_weights_np = np.polynomial.legendre.leggauss(n_r_points)
            u = jnp.asarray((u_np + 1) / 2, dtype=dtype).reshape((1, n_r_points))
            u_weights = jnp.asarray(u_weights_np, dtype).reshape((1, n_r_points))

            segment_index = jnp.arange(n_segments, dtype=jnp.int32).reshape((n_segments, 1))

            grid = {'u': u, 'theta': jnp.array([0]), 'weights': u_weights, 'segment_index': segment_index}
        else:
            u_np, u_weights_np = np.polynomial.legendre.leggauss(n_r_points)
            u = jnp.asarray((u_np + 1) / 2, dtype=dtype).reshape((1, n_r_points, 1))
            u_weights = jnp.asarray(u_weights_np, dtype).reshape((1, n_r_points, 1))

            mu_np, mu_weights_np = np.polynomial.legendre.leggauss(n_theta_points)

            mu = jnp.asarray(mu_np, dtype=dtype).reshape((1, 1, n_theta_points))
            theta = jnp.arccos(mu)

            theta_weights = jnp.asarray(mu_weights_np, dtype=dtype).reshape((1, 1, n_theta_points))

            segment_index = jnp.arange(n_segments, dtype=jnp.int32).reshape((n_segments, 1, 1))

            grid = {'u': u, 'theta': theta, 'weights': u_weights * theta_weights,
                    'segment_index': segment_index}

        self._quadrature_grid_cache[key] = grid
        return grid

    def _get_moment_value_and_gradient_fcn(self,
                                           moment: MomentExpr,
                                           r_cut: float = 50.0,
                                           n_r_points: int = 1024,
                                           n_theta_points: int = 128,
                                           dtype=jnp.float64):

        moment_key = moment if isinstance(moment, str) else id(moment)
        key = (moment_key, r_cut, n_r_points, n_theta_points, str(dtype))

        if key in self._moment_vg_cache:
            return self._moment_vg_cache[key]

        moment_expr = compile_moment_expr(moment) if isinstance(moment, str) else moment
        sig = inspect.signature(moment_expr)
        if unknown := set(sig.parameters) - {'r', 'theta'}:
            raise ValueError(f'Moment function has unknown parameters {unknown}, must only contain r and theta.')
        elif 'r' not in sig.parameters:
            raise ValueError(f'Moment function must contain r in the parameters. Only contains {set(sig.parameters)}.')

        moment_uses_theta = True if 'theta' in sig.parameters else False
        angular = moment_uses_theta or self._fcn_type == 'deformed'

        if self._use_radial_segments:
            res = self._get_segmented_moment_value_and_gradient_fcn(moment_expr, angular, moment_uses_theta, r_cut,
                                                                    n_r_points, n_theta_points, dtype)
        else:
            res = self._get_fixed_moment_value_and_gradient_fcn(moment_expr, angular, moment_uses_theta, r_cut,
                                                                n_r_points, n_theta_points, dtype)
        # self._moment_vg_cache[key] = res
        return res

    def _get_segmented_moment_value_and_gradient_fcn(self, moment_expr: Callable, angular: bool,
                                                     moment_uses_theta: bool,
                                                     r_cut: float, n_r_points: int, n_theta_points: int, dtype):

        grid = self._get_segment_reference_grid(angular, n_r_points, n_theta_points, dtype)

        u = grid['u']
        theta = grid['theta']
        ref_weights = grid['weights']
        segment_index = grid['segment_index']

        def core(pars):

            lower, upper = self.radial_segments(theta.flatten(), pars, r_cut, dtype)

            lower = _reshape_segment_bound(lower, angular, dtype)
            upper = _reshape_segment_bound(upper, angular, dtype)

            dr_du = upper - lower
            r = lower + u * dr_du

            weights = ref_weights * dr_du * r * r
            density = self.density_on_segments(r, theta, pars, segment_index)

            mass_weights = weights * density
            norm = jnp.sum(mass_weights)

            if moment_uses_theta:
                moment_values = moment_expr(r=r, theta=theta)
            else:
                moment_values = moment_expr(r=r)

            numerator = jnp.sum(moment_values * mass_weights)
            return numerator / norm

        return jax.jit(jax.value_and_grad(core))

    def _get_fixed_moment_value_and_gradient_fcn(self, moment_expr: Callable, angular: bool, moment_uses_theta: bool,
                                                 r_cut: float, n_r_points: int, n_theta_points: int, dtype):

        key = (angular, moment_uses_theta, r_cut, n_r_points, n_theta_points, dtype)
        if key in self._fixed_vg_cache:
            r, theta, core = self._fixed_vg_cache[key]
        else:
            grid = self._get_segment_reference_grid(angular, n_r_points, n_theta_points, dtype)

            u = grid['u']
            theta = grid['theta']
            ref_weights = grid['weights']
            segment_index = grid['segment_index']

            lower, upper = self.radial_segments(theta.flatten(), None, r_cut, dtype)

            lower = _reshape_segment_bound(lower, angular, dtype)
            upper = _reshape_segment_bound(upper, angular, dtype)

            dr_du = upper - lower
            r = lower + u * dr_du
            weights = ref_weights * dr_du * r * r

            @jax.jit
            def core(pars, moment_vals):
                density = self.density_on_segments(r, theta, pars, segment_index)

                mass_weights = weights * density
                norm = jnp.sum(mass_weights)

                numerator = jnp.sum(moment_vals * mass_weights)
                return numerator / norm

            self._fixed_vg_cache[key] = (r, theta, core)

        if moment_uses_theta:
            moment_values = moment_expr(r=r, theta=theta)
        else:
            moment_values = moment_expr(r=r)

        return jax.tree_util.Partial(jax.value_and_grad(core, argnums=[0]), moment_vals=moment_values)

    def calculate_moment_gradient(self, moment, pars: dict[str, float], r_cut: float = 50.0, n_r_points: int = 1024,
                                  n_theta_points: int = 128):
        v, g = self._get_moment_value_and_gradient_fcn(moment, r_cut=r_cut, n_r_points=n_r_points,
                                                       n_theta_points=n_theta_points)(pars)
        if isinstance(g, tuple) and len(g) == 1:
            g = g[0]

        return v, g


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
        return 1 / (1 + jnp.exp(4 * jnp.log(3) * (r - c * (1 + beta2 * self.sph_harm_20(theta))) / t))

    def sph_harm_20(self, theta):
        return 1 / 4 * jnp.sqrt(5 / jnp.pi) * (3 * jnp.cos(theta) ** 2 - 1)


class ThreeParameterFermiC(ChargeDistribution):
    def __init__(self, name: str = '3pfc', parameter_names: tuple[str, ...] = ('c', 't', 'beta2'),
                 test_r: npt.NDArray | None = None):
        self._fcn_type = 'deformed'
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def __call__(self, r: Array, theta: Array, pars: dict[str, float]) -> Array:
        c, t, beta2 = pars['c'], pars['t'], pars['beta2']
        return 1 / (1 + jnp.exp(4 * jnp.log(3) * (r * (1 + beta2 * (3 * jnp.cos(theta) ** 2 - 1) / 2) - c) / t))


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
        self._use_radial_segments = True
        self._n_radial_segments = 1
        test_r = test_r if test_r is not None else np.array([0, 5, 10])
        super().__init__(name, parameter_names, test_r)

    def radial_segments(self, theta, pars, r_cut: float, dtype):
        r_cutoff = jnp.asarray(pars["R"], dtype=dtype)

        lower = jnp.stack([0.0])
        upper = jnp.stack([r_cutoff])

        return lower, upper

    def density_on_segments(self, r, theta, pars: dict[str, float], segment_index):
        return jnp.where(segment_index == 0.0, 1.0, 0.0) * jnp.ones_like(r)

    def __call__(self, r: Array, _, pars: dict[str, float]) -> Array:
        r_cutoff = pars['R']
        return -jnp.tanh((r - r_cutoff) * 100) + 1


cds = {cd.name: cd for cd in [TwoParameterFermi(), ThreeParameterFermi(), ParabolicFermi(), Uniform()]}

# if __name__ == '__main__':
#     num_runs = 10000
#
#     threepf = ThreeParameterFermi()
#     test_par = {'c': 6.212, 't': 2.18, 'beta2': 0.324}
#     threepf_moment, threepf_gradient = threepf.calculate_moment_gradient('r**2', pars=test_par)
#     print(np.sqrt(threepf_moment), threepf_gradient)
#
#     uniform = Uniform()
#     print(uniform.calculate_moment_gradient('r**2.2*exp(-0.15*r)', pars={'R': 5.0}))
#     exit()
#
#     uniform = Uniform()
#     print(uniform.calculate_moment_gradient('r**2', pars={'R': 5.0}, n_r_points=2048))
#
#     curr = time.time()
#     for i, k in enumerate(jnp.linspace(1.9, 2.2, 100)):
#         for j, alpha in enumerate(jnp.linspace(0.05, 0.1, 100)):
#             test, _ = twopf.calculate_moment_gradient(f'r**{k}*exp(-{alpha}*r)', pars=test_par)
#         jax.block_until_ready(test)
#         print(i * 100,  time.time() - curr)
#
#     start = time.time()
#     vg_func = jax.block_until_ready(twopf._get_moment_value_and_gradient_fcn('r**2'))
#     print(f'{time.time() - start: 0.3g} s to compile the moment function.')
#
#     start = time.time()
#     for _ in range(num_runs):
#         value, gradient = vg_func(test_par)
#     jax.block_until_ready(gradient)
#     post_twopf = time.time()
#     print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
#           f' or {time_taken / num_runs: 0.3g} s per integral.')
#     print(value, gradient)
#
#     print('Same but more trapezoids.')
#     start = time.time()
#     vg_func = jax.block_until_ready(twopf._get_moment_value_and_gradient_fcn('r**2', n_r_points=4096))
#     print(f'{time.time() - start: 0.3g} s to compile the moment function.')
#
#     num_runs = 10000
#     start = time.time()
#     for _ in range(num_runs):
#         value, gradient = vg_func(test_par)
#     jax.block_until_ready(gradient)
#     post_twopf = time.time()
#     print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
#           f' or {time_taken / num_runs: 0.3g} s per integral.')
#     print(value, gradient)
#
#     print('\nAnd now for the deformed case:')
#     threepf = ThreeParameterFermi()
#     test_par = {'c': 5.0, 't': 2.3, 'beta2': 0.3}
#
#     start = time.time()
#     vg_func = jax.block_until_ready(threepf._get_moment_value_and_gradient_fcn('r**2'))
#     print(f'{time.time() - start: 0.3g} s to compile the moment function.')
#
#     start = time.time()
#     for _ in range(num_runs):
#         value, gradient = vg_func(test_par)
#     jax.block_until_ready(gradient)
#     post_twopf = time.time()
#     print(f'Took {(time_taken := (post_twopf - start)): 0.3g} s for {num_runs} integrals,'
#           f' or {time_taken / num_runs: 0.3g} s per integral.')
#     print(value, gradient)
