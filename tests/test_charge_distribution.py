from functools import partial

import numpy as np
import pytest
from numdifftools import Gradient
from scipy.integrate import dblquad, quad

from radbase.charge_distribution import (ChargeDistribution, cds,
                                         compile_moment_expr)

rtol_value = 1e-5
rtol_grad = 1e-3
atol_grad = 1e-6


def compare_to_scipy(cd, moment, angular, pars):

    moment_func = compile_moment_expr(moment)

    def scipy_moment(par_values, par_names):

        params = {par_name: par_val for par_name, par_val in zip(par_names, par_values)}

        def blank(*_):
            return 1.0

        def scipy_integrate(moment_fn):
            if not angular:
                value, err, *_ = quad(lambda r: moment_fn(r) * r ** 2 * cd(r, None, params), 0, 50, epsabs=1e-6)
            else:
                value, err, *_ = dblquad(
                    lambda r, theta: moment_fn(r) * r ** 2 * np.sin(theta) * cd(r, theta, params),
                    0, np.pi, 0, 50, epsrel=1e-8)
            return value

        return scipy_integrate(moment_func) / scipy_integrate(blank)

    par_n, par_v = list(pars.keys()), list(pars.values())
    scipy_value = scipy_moment(par_v, par_n)
    scipy_gradient = Gradient(partial(scipy_moment, par_names=par_n))(par_v)
    scipy_gradient = {par_n: grad_v for par_n, grad_v in zip(par_n, scipy_gradient.reshape((len(par_n))))}

    cd_value, cd_grad = cd.calculate_moment_gradient(moment, pars)

    assert np.allclose(float(cd_value), scipy_value, rtol=rtol_value)
    for par_name, calc_grad in cd_grad.items():
        assert np.allclose(float(calc_grad), scipy_gradient[par_name], rtol=rtol_grad, atol=atol_grad)


def compare_to_ref(cd: ChargeDistribution, moment: str, pars: dict[str, float], expected: tuple[float, dict[str, float]]):
    val, grad = cd.calculate_moment_gradient(moment, pars)
    print(val, grad)
    assert np.allclose(float(val), expected[0], rtol=rtol_value)
    for par_name, calc_grad in grad.items():
        assert np.allclose(float(calc_grad), expected[1][par_name], rtol=rtol_grad, atol=atol_grad)


def test_uniform():
    uniform = cds['uniform']

    uniform_test_sets = [
        ({'R': 4.0}, 'r**2', (9.6, {'R': 4.8})),
        ({'R': 4.0}, 'r**2.2*exp(-0.05*r)', (10.30304, {'R': 5.23658})),
    ]

    for par, moment, expect in uniform_test_sets:
        compare_to_ref(uniform, moment, par, expect)
        # compare_to_scipy(uniform, moment, False, par)


def test_two_parameter_fermi():
    twopf = cds['2pf']

    twopf_test_sets = [
        ({'c': 5.0, 't': 2.3}, 'r**2', (18.78507, {'c': 6.00014, 't': 3.29133})),
        ({'c': 5.0, 't': 2.3}, 'r**2.2*exp(-0.05*r)', (20.07862, {'c': 6.30665, 't': 3.37103})),
    ]

    for par, moment, expect in twopf_test_sets:
        compare_to_ref(twopf, moment, par, expect)
        compare_to_scipy(twopf, moment, False, par)


def test_three_parameter_fermi():
    threepf = cds['3pf']

    threepf_test_sets = [
        ({'c': 5.0, 't': 2.3, 'beta2': 0.0}, 'r**2', (18.78507, {'c': 6.00014, 't': 3.29133, 'beta2': 0})),
        ({'c': 5.0, 't': 2.3, 'beta2': 0.3}, 'r**2', (19.53567, {'c': 6.31643427204, 't': 3.2561618426, 'beta2': 5.11266584884})),
        ({'c': 5.0, 't': 2.3, 'beta2': 0.3}, 'r**2.2*exp(-0.05*r)', (20.8507667437, {'c': 6.60662568434, 't': 3.32121692178, 'beta2': 5.24610576868})),
    ]

    for par, moment, expect in threepf_test_sets:
        compare_to_ref(threepf, moment, par, expect)
        # compare_to_scipy(threepf, moment, True, par)
