from abc import ABC

from .radius_analyzer import RadiiInformation


class ConversionStep(ABC):
    pass


class ConversionStrategy:
    pass


class ConversionRecord:

    def __init__(self):
        pass


class Measurement(ABC):

    def __init__(self):
        pass


class Converter:

    def __init__(self, conversion_strategy: ConversionStrategy):
        self.strategy = conversion_strategy

    def convert(self, measurements: list[Measurement], conversion_steps) -> RadiiInformation:  # or something like RI
        pass
