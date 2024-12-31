from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import product
from typing import Optional, Collection
from lmfit import Parameters, Minimizer
from dataclasses import dataclass
from uncertainties import ufloat

import copy
import re
import numpy as np


class Nuclide:
    z: int
    a: int

    def __init__(self, nucstr):
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

    def __hash__(self):
        return self.z * 1000 + self.a

    def __str__(self):
        return f'R{self.z:03}{self.a:03}'


class Term:
    """
    Returns some function of the input radii, something like 'R001002' or 'R012014-R012012'
    """
    termstr: str
    termtype: str

    def __init__(self, termstr: str):
        self.termstr = termstr
        self.make_termtype()

    def eval(self, params):
        return params.eval(self.termstr)

    def make_termtype(self):
        if re.fullmatch(r'R\d{6}', self.termstr):
            self.termtype = 'absolute'
        elif re.fullmatch(r'R\d{6}-R\d{6}', self.termstr):
            self.termtype = 'linear relative'
        elif re.fullmatch(r'R\d{6}\*\*2-R\d{6}\*\*2', self.termstr):
            self.termtype = 'squared relative'
        else:
            raise NotImplementedError(f'Term string of {self.termstr} does not match any patterns.')

    def get_nuclides(self):
        return [Nuclide(pstr) for pstr in re.findall(r'R\d{6}', self.termstr)]

    def __repr__(self):
        print(f'{self.termtype} Term: {self.termstr}')


class Citation:
    # stores information about papers
    pass


class RadiusData:
    """
    Describes a probability distribution over radii, along with correlations with any other data
    """
    term: Term
    value: float
    unc: float
    data_id: int
    correlations: defaultdict[int, dict]
    nuclides: list[Nuclide]

    def __init__(self, termstr, value, unc, data_id, **kwargs):
        self.term = Term(termstr)
        self.value = value
        self.unc = unc
        self.data_id = data_id
        self.correlations = defaultdict(dict)
        self.nuclides = self.term.get_nuclides()

    def residual(self, params, normalize=True):
        resid = self.term.eval(params=params) - self.value
        if normalize:
            return resid / self.unc
        return resid

    def copy(self):
        return copy.deepcopy(self)


class Measurement(RadiusData):
    references = list[Citation]

    def __init__(self, references, **kwargs):
        super().__init__(**kwargs)
        self.references = references

    @abstractmethod
    def calc_radius(self):
        pass


class Projection(RadiusData):
    connected_measurements = list[Measurement]

    def __init__(self, connected_measurements, **kwargs):
        super().__init__(**kwargs)
        self.connected_measurements = connected_measurements


class RadiiInformation(dict):
    """
    A dict of RadiusData objects indexed by their data_id, and can calculate residuals/chisq from these data points.
    """

    def __init__(self):
        super().__init__()
        self.MAX_MEASUREMENTS = 10000

    def add(self, term, value, unc, data_id=None, dtype=RadiusData, **kwargs):
        """
        Add a RadiusData

        Parameters
        ----------
        term: str
            quantity described by data; examples are:
            'R001002', the radius of H-2 (Deuterium)
            'R006014-R006012', the difference in radius between C-14 and C-12
            'R006014-R006012', the difference in squared radius between C-14 and C-12
        value: float
            mean value of the `term` quantity
        unc: float
            uncertainty of the `term` quantity
        data_id: int | None
            identifier of the data object. If not given, the smallest integer not yet taken will be chosen automatically
        dtype: type[RadiusData]
            the class of this data; for example, a measurement or a projection
        """
        if data_id in self:
            raise ValueError(f'data_id {data_id} is already in the radii information')

        if data_id is None:
            i = 0
            while i in self:
                if i >= self.MAX_MEASUREMENTS:
                    raise Warning(f'Radii Information has exceeded maximum measurements')
                i += 1
            data_id = i

        new_data = dtype(termstr=term, value=value, unc=unc, data_id=data_id, **kwargs)
        self[data_id] = new_data

    def _get_data_sorted(self):
        return sorted(self.values(), key=lambda d: d.data_id)

    def residuals(self, params, normalize=True):
        return np.array([data.residual(params, normalize) for data in self._get_data_sorted()])

    def get_covariance_matrix(self) -> np.ndarray:
        """
        Generates covariance matrix from correlations of RadiusData objects.
        cov[i][j] is the covariance of ith and jth RadiusData when sorted in ascending data_id order

        Returns
        -------
        np.ndarray:
            Covariance matrix

        """
        data = self._get_data_sorted()
        unc = np.diag([d.unc for d in data])
        correlation = self.get_correlation_matrix()
        return unc.T @ correlation @ unc

    def get_correlation_matrix(self) -> np.ndarray:
        data = self._get_data_sorted()
        ids = [d.data_id for d in data]
        correlation = np.eye(len(data))
        for d in data:
            for correl_id in d.correlations:
                if correl_id not in self:
                    continue
                correlation[ids.index(d.data_id)][ids.index(correl_id)] = d.correlations[correl_id]
        return correlation

    def join(self, other: Collection):
        new_rinfo = self.copy()
        for ri in other:
            new_rinfo = new_rinfo | ri
        return new_rinfo

    def correlate(self, id1: int, id2: int, correlation: float):
        if correlation > 1 or correlation < -1:
            raise ValueError('Correlation must be between -1 and 1')
        if id1 == id2:
            raise ValueError('Cannot correlate data with itself')
        self[id1].correlations[id2] = correlation
        self[id2].correlations[id1] = correlation

    @staticmethod
    def split(datalist: list[RadiusData]):
        new_radinfo = RadiiInformation()
        for data in datalist:
            new_radinfo[data.data_id] = data.copy()

        return new_radinfo

    @staticmethod
    def from_params(params: Parameters, data_kwargs: Optional[dict[str, dict]] = None,
                    dtype: type[RadiusData] = RadiusData):

        rad_info = RadiiInformation()
        for par in params:
            rad_info.add(par, params[par].value, params[par].unc, dtype=dtype, **data_kwargs[par])
        return rad_info

    @property
    def nuclides(self) -> list[Nuclide]:
        return list(set().union(*[data.nuclides for data in self.values()]))


class RadiusAnalyzer:

    def __init__(self):
        self.default_mini_kwargs = {'method': 'lstsq', 'scale_covar': False, 'calc_covar': True}

    def optimize_radii(self, rad_info: RadiiInformation, params: Optional[Parameters] = None,
                       mini_kwargs: Optional[dict] = None) -> RadiiInformation:
        # Same logic defined in solve.py right now
        if params is None:
            params = self.guess_params(rad_info)

        if mini_kwargs is None:
            mini_kwargs = {}
        mini_kwargs = self.default_mini_kwargs | mini_kwargs  # override default mini_kwargs

        mini = Minimizer(rad_info.residuals, params=params, **mini_kwargs)
        mini_res = mini.minimize()
        mini_res.params = getattr(mini_res, 'params')  # purely for type hinting purposes

        new_rad_info = RadiiInformation.from_params(params=mini_res.params)

        return new_rad_info

    def combine_redundant(self, rad_info: RadiiInformation):
        # combines any redundant measurements by taking a weighted average of any data points with the same Term.
        # The combined measurements will also be correlated with each other, something missing from optimize_radii
        # as of 12/29/2024
        terms = set([data.term for data in rad_info])  # collect unique terms

    @staticmethod
    def guess_params(rad_info: RadiiInformation) -> Parameters:
        params = Parameters()
        for nuc in rad_info.nuclides:
            params.add(str(nuc), nuc.guess_radius())
        return params


class DataGrouper:

    @staticmethod
    def group(rad_info: RadiiInformation) -> list[RadiiInformation]:

        def build_graph() -> tuple[dict, dict]:
            edges = defaultdict(set)
            nodes = defaultdict(set)

            for data in rad_info.values():
                nuclides = data.nuclides
                for nuc1, nuc2 in product(nuclides, nuclides):
                    if nuc1 == nuc2:
                        nodes[nuc1].add(data)
                    else:
                        edges[nuc1].add(nuc2)
                        edges[nuc2].add(nuc1)
            return edges, nodes

        def graph_bfs(edges, nodes) -> list[RadiiInformation]:
            groups = []
            seen = set()
            for nuc in nodes:
                if nuc in seen:
                    continue
                group = set()
                stack = [nuc]
                while stack:  # Floodfill / BFS
                    curr_nuc = stack.pop()
                    if curr_nuc in seen:
                        continue
                    seen.add(curr_nuc)

                    group = group.union(nodes[curr_nuc])  # add measurements that touch this nuclide
                    stack += [next_nuc for next_nuc in edges[curr_nuc]]

                groups.append(rad_info.split(list(group)))

            return groups

        return graph_bfs(*build_graph())

    def group_by_lab(self, rad_info: RadiiInformation) -> list[RadiiInformation]:
        pass

    @staticmethod
    def group_by_isotope(rad_info: RadiiInformation) -> list[RadiiInformation]:
        z_groups = defaultdict(list)
        for data in rad_info.values():
            atomic_numbers = [nuc.z for nuc in data.nuclides]
            if len(set(atomic_numbers)) > 1:  # this measurement extends between isotopes
                continue
            z_groups[atomic_numbers[0]].append(data)

        z_groups = sorted([(atomic_number, rad_info.split(data)) for atomic_number, data in z_groups.items()])
        return [group[1] for group in z_groups]


@dataclass
class RadiusDataRequest:
    # contains what consumer wants to know as well as filters on received data
    nuclides: list[str] | None


@dataclass
class RadiusDataResponse:
    # response sent back, can contain measurement, analysis, and bibliographic information
    # instances will be copied versions of original variables, so mutate away!
    measurement_info: RadiiInformation
    optimized_radii: RadiiInformation


class RadiusDataRequester(ABC):

    @abstractmethod
    def respond(self, request: RadiusDataRequest) -> RadiusDataResponse:
        pass


class RadiusDataGateway(ABC):

    @abstractmethod
    def load_data(self, **kwargs) -> RadiiInformation:
        pass

    @abstractmethod
    def load_analysis(self, **kwargs) -> RadiiInformation:
        pass


class RadiusDataGenerator(RadiusDataRequester):

    def __init__(self, data_loader: RadiusDataGateway):
        self.data_loader = data_loader
        self.measurement_information = self.data_loader.load_data()
        self.analysis_information = self.data_loader.load_analysis()

    @staticmethod
    def filter(request: RadiusDataRequest, rad_info: RadiiInformation) -> RadiiInformation:
        return rad_info.split([data for data in rad_info.values() if any(nuc in request.nuclides for nuc in data)])

    def parse_request(self, request: RadiusDataRequest):
        meas_info = self.filter(request, self.measurement_information)
        analysis_info = self.filter(request, self.analysis_information)

        return meas_info, analysis_info

    @staticmethod
    def convert_to_response(meas_info, analysis_info):
        return RadiusDataResponse(meas_info, analysis_info)

    def respond(self, request: RadiusDataRequest):
        return self.convert_to_response(*self.parse_request(request))
