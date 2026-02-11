"""
Measurement -> ConvertedQuantity -> ... -> Radius

TODO: Consider better way of defining the type of a quantity. For example, each quantity will have a value and an
uncertainty, and then some kind of identifier (just a string?) for what those numbers go with. For now, I will use the
naive approach of using strings, although this will definitely come back to haunt me.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Callable, Literal, Tuple

from uncertainties import UFloat, ufloat

# from .radius_analyzer import Nuclide

# from .radius_analyzer import RadiiInformation


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


class Transition:

    def __init__(self):
        raise NotImplementedError()


class QuantityType(str):
    pass


@dataclass
class Quantity:
    uvar: UFloat
    value: float = field(init=False)
    unc: float = field(init=False)
    unit: str
    quantity: QuantityType

    def __post_init__(self):
        self.value = self.uvar.n
        self.unc = self.uvar.s


class MuonicTransitionEnergyMeasurement:
    reference: str
    transition: Transition
    nuclide: Nuclide
    uvar: UFloat
    unit: Literal['eV']


class OpticalIsotopeShiftMeasurement:
    reference: str
    transition: Transition
    nuclide_heavier: Nuclide
    nuclide_lighter: Nuclide
    uvar: UFloat
    unit: Literal['THz']


@dataclass
class ConversionHistory:

    """Stores one step of how a converted quantity was obtained. Contains input quantites and ConverstionStep(s) followed
    to reach the ConvertedQuantity. By following the chain of ConversionHistories, the full conversion process can
    be reconstructed."""

    conversion_text: str
    input_quantities: Tuple[Quantity]
    output_quantity: Quantity
    conversion_func: Callable | None


@dataclass
class Measurement(Quantity):
    ref: str | None


@dataclass
class ConvertedQuantity(Quantity):
    history: ConversionHistory


class ConversionStep(ABC):
    """
    Wrapper around a conversion function that indicates input data types, creates the conversion history for the output,
    and returns the output type.
    """

    input_quantities: tuple[QuantityType]
    output_quantity: QuantityType
    conversion_function: Callable[[*Tuple[ufloat]], ufloat]  # Variable number of ufloats as input, returns ufloat as output

    def __init__(self, input_quantities: tuple[QuantityType], output_quantity: QuantityType,
                 conversion_function: Callable[[*Tuple[ufloat]], ufloat]):
        self.input_quantities = input_quantities
        self.output_quantity = output_quantity
        self.conversion_function = conversion_function

    @abstractmethod
    def convert_quantity(self, *inputs) -> ConvertedQuantity:
        pass

    @abstractmethod
    def create_history(self, *inputs) -> ConversionHistory:
        pass


class ConversionStrategy:
    pass


class Converter:

    def __init__(self, conversion_strategy: ConversionStrategy):
        self.strategy = conversion_strategy

    def convert(self, measurements: list[Measurement], conversion_steps: list[Quantity]):  # or something like RI
        pass
