from abc import ABC, abstractmethod


class Measurement(ABC):

    @abstractmethod
    def calc_radius(self, information):
        pass

    @abstractmethod
    def calc_seltzer(self):
        pass


class MuonicMeasurement(Measurement):
    pass


class ElectronicMeasurement(Measurement):
    pass


class KalphaMeasurement(Measurement):
    pass


class OISMeasurement(Measurement):
    pass


class RadiusConverter:
    measurements = list[Measurement]

    def request_information(self, information):
        """Retrieves desired quantity from measurements; if quantity not found yields a calculated value"""
        pass
