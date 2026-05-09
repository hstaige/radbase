"""
RadQuantity
"""

from functools import partial
from typing import Annotated, Any, Callable, Self, Tuple, TypeAlias

from pint import UnitRegistry
from pint.registry import Quantity as PintQuantity
from sympy import Symbol, lambdify
from sympy.parsing.sympy_parser import parse_expr
from uncertainties import ufloat

UQuantity: TypeAlias = Annotated[
    PintQuantity[Any],
    "pint Quantity with UFloat magnitude.",
]


class Nuclide:
    z: int
    a: int

    def __init__(self, nucstr):
        self.nucstr = nucstr
        nucstr = nucstr[1:]
        self.z = int(nucstr[:3])
        self.a = int(nucstr[3:])
        self.n = self.a - self.z

    def guess_radius(self):
        # semiempirical estimate
        return (0.9071 + 1.1025 / (self.a ** (2 / 3)) + -0.548 / (self.a ** (4 / 3))) * self.a ** (1 / 3)

    def __repr__(self):
        return f'Nuclide (Z={self.z}, A={self.a})'

    def __eq__(self, other):
        return (self.z == other.z) and (self.a == other.a)

    def __lt__(self, other):
        if self.z != other.z:
            return self.z < other.z
        return self.a < other.a

    def __hash__(self):
        return self.z * 1000 + self.a

    def __str__(self):
        return f'R{self.z:03}{self.a:03}'


class ConversionStep:
    conversion_function: Callable[[*Tuple[UQuantity, ...]], UQuantity]
    description: str | None = None

    def __init__(self, conversion_function, description):
        self.conversion_function = conversion_function
        self.description = description

    def __call__(self, **inputs):
        return self.conversion_function(**inputs)

    def __repr__(self):
        return f'ConversionStep: {self.description}'


class ConversionStrategy:
    def value_strategy(self, quantity):
        return quantity.reported_value


conversion_strategy = ConversionStrategy()


class RadQuantity:
    key: str

    reported_value: UQuantity | None = None

    conversion_step: ConversionStep | None = None
    inputs: dict[str, UQuantity | Self] | None  # Variable number of RadQuantities

    metadata: dict | None = None

    def __init__(self, key: str,
                 reported_value: UQuantity | None = None,
                 conversion_step: ConversionStep | Callable[[Self], ConversionStep] | None = None,
                 inputs: dict[str, UQuantity | Self] | None = None,
                 metadata: dict[str, Any] = None,
                 self_inputs: dict[str, str] | None = None):
        self.key = key
        self.reported_value = reported_value
        self.metadata = metadata

        no_inputs = inputs is None and self_inputs is None
        if not (conversion_step is None and no_inputs or conversion_step is not None and not no_inputs):
            raise ValueError(f'{conversion_step=} and {inputs=}, either both or neither must be None.')

        self.inputs = inputs
        if self_inputs:
            self.inputs = self.inputs if self.inputs is not None else {}
            for key, value in self_inputs.items():
                curr = self
                for subvalue in value.split('.'):
                    if subvalue == 'n':
                        curr = q_nominal(curr)
                    else:
                        curr = curr.__getattribute__(subvalue)
                self.inputs[key] = curr

        if conversion_step is not None:
            self.conversion_step = conversion_step if isinstance(conversion_step, ConversionStep) else conversion_step(
                self)
        else:
            self.conversion_step = None

    @property
    def evaluated_value(self) -> UQuantity | None:
        if self.conversion_step is None or self.inputs is None:
            return None
        return self.conversion_step(**self.inputs)

    @property
    def value(self) -> UQuantity:
        return conversion_strategy.value_strategy(self)

    def __repr__(self):
        ret = 'RadQuantity('
        ret += f'key: {self.key=}, '
        ret += f'reported_value: {self.reported_value}, '
        ret += f'evaluated_value: {self.evaluated_value}, '
        ret += ')'
        return ret

    def __hash__(self):
        return self.key

    def __add__(self, other):
        return self.value + getattr(other, 'value', other)

    def __sub__(self, other):
        return self.value - getattr(other, 'value', other)

    def __mul__(self, other):
        return self.value * getattr(other, 'value', other)

    def __truediv__(self, other):
        return self.value / getattr(other, 'value', other)


def make_algebraic_expression(expression: str, description: str = '') -> ConversionStep:
    def algebraic_expression(sympy_expr, **values: dict[str, UQuantity]):
        return sympy_expr(**values)

    expression = parse_expr(expression)
    symbols = tuple(expression.atoms(Symbol))
    conv_function = partial(algebraic_expression, lambdify(symbols, expression, 'numpy'))

    return ConversionStep(conversion_function=conv_function, description=f'{description}: {expression}')


barrett_radius_linear = make_algebraic_expression('rka0 + cz * (energy - energy0)')
barrett_to_radius = make_algebraic_expression('rka / v2')


def make_fixed_conversion(rad_quantity: RadQuantity, description: str = 'Fixed conversion step.') -> ConversionStep:
    """
    Used when a conversion function is opaque and cannot be called with different inputs. Commonly seen when
    doing fits to charge distributions, as you do not know what the output parameters would be for different functions.
    """

    def fixed_function(inputs, output, **values: dict[str, UQuantity]):
        if inputs == values:
            return output
        raise ValueError('Tried to call fixed function with different inputs.')

    conv_function = partial(fixed_function, rad_quantity.inputs, rad_quantity.reported_value)
    return ConversionStep(conversion_function=conv_function, description=description)


def make_cd_fit_conversion(cd_type: str, out_param: str, fit_description: str = '') -> Callable[[RadQuantity], ConversionStep]:
    return partial(make_fixed_conversion, description=f'Parameter {out_param} of {cd_type} fit. {fit_description}')


def q_nominal(q: UQuantity) -> PintQuantity[float]:  # Extracts just the nominal value from a UQuantity
    return q.m.n * q.u


def rq_nominal(rq: RadQuantity) -> PintQuantity[float]:
    return rq.value.m.n * rq.value.u


if __name__ == '__main__':
    ureg = UnitRegistry()

    energy = RadQuantity('test_Uo_500_muonicenergy', ufloat(500, 0.1) * ureg('keV'))

    rka = RadQuantity('test_Uo_500_barrettradius', ufloat(4, 0.1) * ureg('fm'),
                      conversion_step=barrett_radius_linear,
                      self_inputs={'rka0': 'reported_value.n'},
                      inputs={'energy': energy,
                              'cz': ufloat(1, 0.001) * ureg('fm/keV'),
                              'energy0': energy},
                      metadata={'k': 2.3,
                                'alpha': 0.6})

    print(rka.evaluated_value)

    cross_section_q1 = RadQuantity('test_Uo_500_escatt_q1', ufloat(300, 0.1) * ureg('fm^2/steradian'))
    cross_section_q2 = RadQuantity('test_Uo_500_escatt_q2', ufloat(100, 0.1) * ureg('fm^2/steradian'))

    escatt_fit_inputs = (cross_section_q1, cross_section_q2)
    c = RadQuantity('test_Uo_500_fermi_c', ufloat(3, 0.2) * ureg('fm'),
                    conversion_step=make_cd_fit_conversion('2pf', 'c'),
                    inputs={'q1': cross_section_q1,
                            'q2': cross_section_q2})

    t = RadQuantity('test_Uo_500_fermi_c', ufloat(3, 0.2) * ureg('fm'),
                    conversion_step=make_cd_fit_conversion('2pf', 't'),
                    inputs={'q1': cross_section_q1,
                            'q2': cross_section_q2})

    print(c, t, c.conversion_step, t.conversion_step)

    v2 = RadQuantity('test_Uo_500_v2', ufloat(1.2, 0.1),
                     conversion_step=make_fixed_conversion,
                     inputs={'c': c,
                             't': t,
                             'k': rka.metadata['k'],
                             'alpha': rka.metadata['alpha']})

    print(v2.evaluated_value)

    radius = RadQuantity('test_Uo_500_radius', ufloat(3, 0.1) * ureg('fm'),
                         conversion_step=barrett_to_radius,
                         inputs={'rka': rka,
                                 'v2': v2})

    print(radius.evaluated_value)
