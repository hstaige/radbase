import copy
import itertools
import re
from abc import ABC, abstractmethod
from collections import defaultdict, deque
from dataclasses import dataclass
from functools import partial
from itertools import product
from typing import Optional

import numpy as np
import pandas as pd
from lmfit import Minimizer, Parameters
from uncertainties import correlated_values_norm, correlation_matrix, ufloat


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


class Term:
    """
    Returns some function of the input radii, something like 'R001002' or 'R012014-R012012'
    """
    termstr: str
    termtype: str

    def __init__(self, termstr: str):
        self.termstr = termstr
        self.termtype = ''

    @abstractmethod
    def eval(self, params: Parameters) -> ufloat:
        pass

    def asteval(self, params):
        return params.eval(self.termstr)

    def get_nuclides(self):
        return [Nuclide(pstr) for pstr in re.findall(r'R\d{6}', self.termstr)]

    @abstractmethod
    def calc_uvar(self, uvars: dict[str, ufloat]) -> ufloat:
        """Calculates value from dict of uvars"""
        pass

    @abstractmethod
    def adjust_params(self, params):
        pass

    def __repr__(self):
        return f'{self.termtype} Term: {self.termstr}'


def create_term(termstr):
    if re.fullmatch(r'R\d{6}', termstr):
        return AbsoluteTerm(termstr)
    elif re.fullmatch(r'R\d{6}-R\d{6}', termstr):
        return LinearRelativeTerm(termstr)
    elif re.fullmatch(r'R\d{6}\*\*2-R\d{6}\*\*2', termstr):
        return SquaredRelativeTerm(termstr)
    else:
        raise NotImplementedError(f'Term string of {termstr} does not match any patterns.')


class AbsoluteTerm(Term):

    def __init__(self, termstr):
        super().__init__(termstr)
        self.termtype = 'absolute'
        self.rstr = termstr

    def eval(self, params):
        return params[self.rstr].value

    def calc_uvar(self, uvars):
        return uvars[self.rstr]

    def adjust_params(self, params):
        return params


class LinearRelativeTerm(Term):
    """A linear relative term of the form R1-R2"""

    def __init__(self, termstr):
        super().__init__(termstr)
        self.termtype = 'linear_relative'
        self.R1, self.R2 = re.findall(r'R\d{6}', termstr)

    def eval(self, params):
        return params[self.R1].value - params[self.R2].value

    def calc_uvar(self, uvars: dict[str, ufloat]) -> ufloat:
        return uvars[self.R1] - uvars[self.R2]

    def adjust_params(self, params):
        smallest = min(params.keys())
        params[smallest].vary = False
        return params


class SquaredRelativeTerm(Term):

    def __init__(self, termstr):
        super().__init__(termstr)
        self.termtype = 'squared_relative'
        self.R1, self.R2 = re.findall(r'R\d{6}', termstr)

    def eval(self, params):
        return params[self.R1].value ** 2 - params[self.R2].value ** 2

    def calc_uvar(self, uvars: dict[str, ufloat]) -> ufloat:
        return uvars[self.R1] ** 2 - uvars[self.R2] ** 2

    def adjust_params(self, params):
        smallest = min(params.keys())
        params[smallest].vary = False
        return params


class GeneralTerm(Term):
    """
    Still needs serious improvements, this is a placeholder and should not be used in production code.

    A general term to use when one of the existing categories does not fit. Will be slower, as it will use the asteval
    of parameters rather than the params[par].value syntax.
    """
    def __init__(self, termstr):
        super().__init__(termstr)
        self.termtype = 'general'

    def eval(self, params):
        return self.asteval(params)

    def calc_uvar(self, uvars: dict[str, ufloat]) -> ufloat:
        raise NotImplementedError
        return ufloat(0, 0)

    def adjust_params(self, params):
        return params


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

    def __init__(self, term: str | Term, value, unc, data_id, **kwargs):
        if isinstance(term, Term):
            self.term = term
        elif isinstance(term, str):
            self.term = create_term(term)
        else:
            raise TypeError(f'{term} was provided by the type of {type(term)} is not a string or Term')

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

    def __repr__(self):
        return f'Id {self.data_id}: {self.term.termstr} = {self.value} +/- {self.unc}'


class Measurement(RadiusData):

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
        self._uvars = None

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
            data_id = self._get_next_id()

        new_data = dtype(term=term, value=value, unc=unc, data_id=data_id, **kwargs)
        self[data_id] = new_data

    def get_data_sorted(self):
        return sorted(self.values(), key=lambda d: d.data_id)

    def residuals(self, params, normalize=True):
        return np.array([data.residual(params, normalize) for data in self.get_data_sorted()])

    def get_covariance_matrix(self) -> np.ndarray:
        """
        Generates covariance matrix from correlations of RadiusData objects.
        cov[i][j] is the covariance of ith and jth RadiusData when sorted in ascending data_id order

        Returns
        -------
        np.ndarray:
            Covariance matrix

        """
        data = self.get_data_sorted()
        unc = np.diag([d.unc for d in data])
        correlation = self.get_correlation_matrix()
        return unc.T @ correlation @ unc

    def get_correlation_matrix(self) -> np.ndarray:
        data = self.get_data_sorted()
        ids = [d.data_id for d in data]
        correlation = np.eye(len(data))
        for d in data:
            for correl_id in d.correlations:
                if correl_id not in self:
                    continue
                correlation[ids.index(d.data_id)][ids.index(correl_id)] = d.correlations[correl_id]
        return correlation

    def join(self, other):
        """Main challenge is maintaining correlations/correct data ids in this process. For now, simply shift all
        relevant ids UP by a certain amount to avoid overlaps"""
        new_rinfo = copy.deepcopy(self)
        for ri in other:
            ri = ri.copy()
            min_data_id = max(list(new_rinfo.keys()) + [0])
            ri.shift_all_ids(min_data_id + 1)  # move all data ids upwards
            if any(data_id in new_rinfo.keys() for data_id in ri.keys()):
                raise ValueError('Somehow have overlapping keys even after shift!')
            for data in ri.values():
                new_rinfo[data.data_id] = data
        return new_rinfo

    def shift_all_ids(self, amount: int):
        # TODO: Make this less janky. Currently im manually updating the RadiiData,
        #  but this should really be a method in the data object
        if amount == 0:
            return
        overwritten = set()
        for data_id, data in list(self.items()):
            data.data_id += amount
            for other in list(data.correlations.keys()):
                data.correlations[other + amount] = data.correlations[other]
                del data.correlations[other]
            self[data_id + amount] = data
            overwritten.add(data_id + amount)
            if data_id not in overwritten:
                del self[data_id]

    def correlate(self, id1: int, id2: int, correlation: float):
        if id1 == id2:  # don't change correlation matrix for data with itself.
            return
        if correlation > 1 or correlation < -1:
            raise ValueError('Correlation must be between -1 and 1')
        if abs(correlation) < 1e-2:
            correlation = 0
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
        """Constructs a RadiiInformation object from a lmfit.Parameters object."""

        if data_kwargs is None:
            data_kwargs = defaultdict(dict)

        rad_info = RadiiInformation()
        for par in params:
            rad_info.add(par, params[par].value, params[par].stderr, dtype=dtype, **data_kwargs[par])
        return rad_info

    @staticmethod
    def from_terms_and_uvars(terms: list[Term], uvars: list[ufloat], data_kwargs: Optional[dict[str, dict]] = None,
                             dtype: type[RadiusData] = RadiusData):
        """Constructs radii information given terms[i] is described by uvar[i].
        Captures correlations in uvars"""

        if data_kwargs is None:
            data_kwargs = defaultdict(dict)

        rad_info = RadiiInformation()
        for i, (term, uvar) in enumerate(zip(terms, uvars)):
            rad_info.add(term.termstr, uvar.nominal_value, uvar.std_dev, dtype=dtype, **data_kwargs[term.termstr])
        term_to_id = {data.term.termstr: data.data_id for data in rad_info.values()}

        corr_matrix = correlation_matrix(uvars)
        for i, row in enumerate(corr_matrix):
            for j, val in enumerate(row):
                if i == j:
                    continue
                id1, id2 = term_to_id[terms[i].termstr], term_to_id[terms[j].termstr]
                rad_info.correlate(id1, id2, val)

        return rad_info

    def copy(self):
        return copy.deepcopy(self)

    def _get_next_id(self):
        i = 0
        while i in self:
            if i >= self.MAX_MEASUREMENTS:
                raise Warning('Radii Information has exceeded maximum measurements')
            i += 1
        return i

    def create_uvars(self, use_cached=False):
        """Return a list of ufloats that describe the values and covariance matrix of the data. Note that
        uvars correspond to the indexes in ascending order."""
        if use_cached and self._uvars is not None:
            return self._uvars

        val_uncs = [(data.value, data.unc) for data in self.get_data_sorted()]
        corr_matrix = self.get_correlation_matrix()
        self._uvars = correlated_values_norm(val_uncs, corr_matrix)
        return self._uvars

    def evaluate_terms(self, terms: list[str] | list[Term], use_cached=False):
        # if not all(data.term.termtype == 'absolute' for data in self.values()):
        #     raise ValueError('self must only have absolute terms.')

        terms = [create_term(term) if isinstance(term, str) else term for term in terms]
        abs_uvars = {data.term.termstr: abs_uvar for data, abs_uvar in zip(self.get_data_sorted(),
                                                                           self.create_uvars(use_cached=use_cached))}
        term_uvars = [term.calc_uvar(abs_uvars) for term in terms]
        return RadiiInformation.from_terms_and_uvars(terms, term_uvars)

    @property
    def nuclides(self) -> list[Nuclide]:
        return list(set().union(*[data.nuclides for data in self.values()]))


def compare_radii_information(outfile, rad_infos: list[RadiiInformation], prefixes: None | list[str] = None):
    """Dumps term-by-term comparison of rad_infos into an Excel file. Expects each term to show up only once in each
    radii_information object"""

    if prefixes is None:
        prefixes = [str(f'R{i}') for i in range(len(rad_infos))]

    if len(rad_infos) != len(prefixes):
        raise ValueError(f'rad_infos and prefixes must have same length,'
                         f' instead have lengths {len(rad_infos)} and {len(prefixes)}')

    terms = sorted(list(set([data.term.termstr for rad_info in rad_infos for data in rad_info.values()])))
    info = defaultdict(list)
    for term in terms:
        info['Term'].append(term)
        for prefix, rad_info in zip(prefixes, rad_infos):
            match = [data for data in rad_info.values() if data.term.termstr == term]
            if len(match) == 0:
                try:
                    term_data = rad_info.evaluate_terms([term], use_cached=True)[0]
                    info[prefix + '_value'].append(term_data.value)
                    info[prefix + '_unc'].append(term_data.unc)
                except (
                        ValueError,
                        KeyError):  # value error if radinfo not only absolute, key error if radius not found.
                    info[prefix + '_value'].append(np.nan)
                    info[prefix + '_unc'].append(np.nan)
                continue
            # if len(match) > 1:
            #     raise Warning(f'Multiple data for prefix {prefix} have term "{term}"')
            info[prefix + '_value'].append(match[0].value)
            info[prefix + '_unc'].append(match[0].unc)

        if len(rad_infos) == 2:
            info[f'{prefixes[0]}-{prefixes[1]}_value'].append(
                info[prefixes[0] + '_value'][-1] - info[prefixes[1] + '_value'][-1])
            info[f'{prefixes[0]}-{prefixes[1]}_unc'].append(
                info[prefixes[0] + '_unc'][-1] - info[prefixes[1] + '_unc'][-1])

    writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
    format1 = writer.book.add_format({'num_format': '0.0000'})
    df = pd.DataFrame(info)
    df.to_excel(writer, sheet_name='comparison', index=False, na_rep='NaN')

    for column in df:
        column_length = max(df[column].astype(str).map(len).max(), len(column))
        col_idx = df.columns.get_loc(column)
        writer.sheets['comparison'].set_column(col_idx, col_idx, column_length, format1)

    if len(rad_infos) == 2:  # add some conditional formatting to hightlight odd cells
        format_red = writer.book.add_format({'bg_color': '#FFC7CE',
                                             'font_color': '#9C0006'})

        for idx in [len(df.columns) - 2, len(df.columns) - 1]:
            writer.sheets['comparison'].conditional_format(1, idx, len(df), idx, {'type': 'cell',
                                                                                  'criteria': 'not between',
                                                                                  'minimum': '=-MIN(ABS($C2), ABS($E2))',
                                                                                  'maximum': '=MIN(ABS($C2), ABS($E2))',
                                                                                  'format': format_red})

    writer.close()  # (writer.save() was deprecated and removed as of 2023/2024)


class DataGroup:

    def __init__(self, all_data: RadiiInformation, primary_data: Optional[RadiiInformation] = None,
                 secondary_data: Optional[list[RadiiInformation]] = None, name=''):
        self.all_data = all_data
        self.primary_data = primary_data if primary_data else all_data
        self.secondary_data = secondary_data if secondary_data is not None else []
        self.name = name

    def join(self, others):
        new_group = copy.deepcopy(self)
        new_group.all_data = new_group.all_data.join([other.all_data for other in others])
        new_group.primary_data = new_group.primary_data.join([other.primary_data for other in others])
        new_group.secondary_data = list(
            itertools.chain.from_iterable(self.secondary_data + [other.secondary_data for other in others]))
        return new_group

    @staticmethod
    def from_list(data_groups, name=''):
        all_data = RadiiInformation().join([group.all_data for group in data_groups])
        primary_data = RadiiInformation().join([group.primary_data for group in data_groups])
        secondary_data = list(itertools.chain.from_iterable([group.secondary_data for group in data_groups]))
        return DataGroup(all_data, primary_data, secondary_data, name=name)

    @property
    def nuclides(self):
        return self.all_data.nuclides


class RadiusAnalyzer:

    def __init__(self):
        self.default_mini_kwargs = {'scale_covar': False, 'calc_covar': True, 'tol': 1e-5, 'max_nfev': 1e6}

    def optimize_radii(self, data: RadiiInformation | DataGroup, params: Optional[Parameters] = None,
                       method='lbfgsb', mini_kwargs: Optional[dict] = None) -> RadiiInformation:
        """
        Takes in radius data and returns the optimal radius values for each nuclide.

        Note that if data contains only one kind of term (e.g. only linear relative or squared relative),
        the returned RadiiInformation will also contain only that kind of term. This is because if we only know
        'R001002-R001001', 'R001003-R001002', and 'R001003-R001001' we cannot determine individual radii, but we can
        determine differences.

        Parameters
        ----------
        data: RadiiInformation | DataGroup
            Contains radius data to use for optimization
        params: Parameters
            initial values for optimization
        method: str
            The optimization method for lmfit.minimize
        mini_kwargs: dict
            keywords for lmfit.Minimizer

        Returns
        -------
        opt_rad_info: RadiiInformation
        """
        if isinstance(data, RadiiInformation):
            data = DataGroup(all_data=data, primary_data=data)

        if params is None:
            params = self.guess_params(data.all_data)

        params = params.copy()
        for par in list(params.keys()):
            if Nuclide(par) not in data.nuclides:
                del params[par]

        if len(set(data.term.termtype for data in data.all_data.values())) == 1:  # all one type of term
            params = list(data.all_data.values())[0].term.adjust_params(params)

        orig_params = params.copy()

        if mini_kwargs is None:
            mini_kwargs = {}
        mini_kwargs = self.default_mini_kwargs | mini_kwargs  # override default mini_kwargs

        """
        Somewhat annoyingly, lmfit.minimize will not provide uncertainties if len(residuals) <= len(params). To get around this,
        and allow for 'redundant' minimizations, I append a zero to the residuals array, and adjust the shape of the covariance
        matrix to reflect this. This is definitely a hack, but it is easier then checking all edges cases.
        """

        def add_zero_resid(pars, rad_info):
            return np.pad(rad_info.residuals(pars, normalize=False), (0, 1), 'constant')

        def reduce_fcn(resids: np.ndarray, inv_cov: np.ndarray):
            return resids.T @ inv_cov @ resids

        inv_cov_matrix = np.linalg.inv(data.primary_data.get_covariance_matrix())
        inv_cov_matrix = np.pad(inv_cov_matrix, (0, 1), 'constant')
        # noinspection PyTypeChecker
        mini_kwargs['reduce_fcn'] = partial(reduce_fcn, inv_cov=inv_cov_matrix)

        # Fits the primary data, keeping any nuclides that don't affect primary data constant
        for p in params:
            if Nuclide(p) not in data.primary_data.nuclides:
                params[p].vary = False

        mini = Minimizer(add_zero_resid, params=params, fcn_args=(data.primary_data,), **mini_kwargs)
        mini_res = mini.minimize(method=method)
        mini_res.params = getattr(mini_res, 'params')  # purely for type hinting purposes
        mini_res.covar = getattr(mini_res, 'covar')
        params = mini_res.params

        for p in params:
            params[p].vary = False

        for secondary_info in data.secondary_data:
            for p in params:
                if (Nuclide(p) in secondary_info.nuclides) and (Nuclide(p) not in data.primary_data.nuclides):
                    params[p].vary = True

            inv_cov_matrix = np.linalg.inv(secondary_info.get_covariance_matrix())
            inv_cov_matrix = np.pad(inv_cov_matrix, (0, 1), 'constant')
            # noinspection PyTypeChecker
            mini_kwargs['reduce_fcn'] = partial(reduce_fcn, inv_cov=inv_cov_matrix)

            mini = Minimizer(add_zero_resid, params=params, fcn_args=(secondary_info,), **mini_kwargs)
            mini_res = mini.minimize(method=method)
            params = mini_res.params

            for p in params:
                params[p].vary = False

        for p in params:
            params[p].vary = orig_params[p].vary

        inv_cov_matrix = np.linalg.inv(data.all_data.get_covariance_matrix())
        inv_cov_matrix = np.pad(inv_cov_matrix, (0, 1), 'constant')
        # noinspection PyTypeChecker
        mini_kwargs['reduce_fcn'] = partial(reduce_fcn, inv_cov=inv_cov_matrix)

        mini = Minimizer(add_zero_resid, params=params, fcn_args=(data.all_data,), **mini_kwargs)
        mini_res = mini.minimize(method=method)
        terms = [create_term(par) for par in mini_res.params.keys()]
        uvars = mini_res.params.create_uvars(covar=mini_res.covar)

        if len(set(dat.term.termtype for dat in data.all_data.values())) == 1:  # all one type of term
            terms = list(set([dat.term for dat in data.all_data.values()]))
            uvars = {term: term.calc_uvar(uvars) for term in terms}

        uvars = list(uvars.values())
        new_rad_info = RadiiInformation.from_terms_and_uvars(terms, uvars)

        return new_rad_info

    def adjust_uncertainties(self, rad_info):
        optimized_rad_info = RadiiInformation().join([self.optimize_radii(group) for group in
                                                      DataGrouper.group_by_term(rad_info)])
        return self.adjust_uncertainties_by_term(rad_info, optimized_rad_info)

    @staticmethod
    def adjust_uncertainties_by_term(original_rad_info: RadiiInformation,
                                     optimized_rad_info: RadiiInformation) -> RadiiInformation:
        """

        Parameters
        ----------
        original_rad_info: Original RadiiInformation with all data
        optimized_rad_info: An optimized RadiiInformation, potenntially containing only a subset of the data

        Returns
        -------
        """
        term_data = {data.term.termstr: data for data in optimized_rad_info.values()}
        for term, opt_data in term_data.items():

            match = [data for data in original_rad_info.values() if data.term.termstr == term]
            if len(match) <= 1:
                continue
            values, uncs = list(zip(*[(data.value, data.unc) for data in match]))
            values, uncs = np.array(values), np.array(uncs)
            red_chi_sq = np.sum(np.square((values - opt_data.value) / uncs)) / (len(match) - 1)
            multiplier = max(1, np.sqrt(0.5 * (1 + red_chi_sq)))
            for data in match:
                original_rad_info[data.data_id].unc *= multiplier
        return original_rad_info

    @staticmethod
    def guess_params(rad_info: RadiiInformation) -> Parameters:
        params = Parameters()
        for nuc in rad_info.nuclides:
            params.add(str(nuc), nuc.guess_radius())
        return params

    def adjust_from_previous_optimization(self, original_rad_info: RadiiInformation,
                                          optimized_rad_info: RadiiInformation) -> Parameters:
        """If additional nuclides are present in original rad info that are not in optimized rad info"""

        init_params = self.guess_params(original_rad_info)
        for par in init_params.values():
            par.vary = False
        for data in optimized_rad_info.values():
            if data.term.termtype == 'absolute':
                init_params[data.term.termstr].value = data.value

        for group in DataGrouper().group_by_element(original_rad_info):
            if all(nuc in optimized_rad_info.nuclides for nuc in group.nuclides):
                continue

            for nuc in group.nuclides:
                init_params[nuc.nucstr].vary = nuc not in optimized_rad_info.nuclides

            new_rad_info = self.optimize_radii(group, params=init_params)
            for data in new_rad_info.values():
                if data.term.termstr == 'absolute':
                    init_params[data.term.termstr] = data.value

            for nuc in group.nuclides:
                init_params[nuc.nucstr].vary = False

        for par in init_params.values():
            par.vary = True

        return init_params

    @staticmethod
    def calculate_f_factors(rad_info: RadiiInformation,
                            opt_rad_info: RadiiInformation) -> tuple[dict[int, ufloat], dict[int, ufloat]]:
        """Calculates the 'f' correction factor for optical isotope shift measurements

        Note for future Hunter: this function uses the data.measurement_type attribute I hacked in in database.py
        should definitely be changed once the conversion has been implemented.
        """
        not_ois_nuclides = set()
        for data in rad_info.values():
            if getattr(data, 'measurement_type', 'ois') != 'ois':
                not_ois_nuclides.update(set(data.nuclides))

        f_factors = {}
        delta_factors = {}
        opt_uvars = opt_rad_info.create_uvars()
        opt_uvars = {data.term.termstr: uvar for data, uvar in zip(opt_rad_info.get_data_sorted(), opt_uvars)}
        for group in DataGrouper.group_by_element(rad_info):
            z = group.nuclides[0].z
            group = group.all_data
            fs, ds = [], []
            for data in group.values():
                try:
                    if getattr(data, 'measurement_type', '') == 'ois' and all(nuc in not_ois_nuclides for nuc in data.nuclides):
                        optimal = data.term.calc_uvar(opt_uvars)
                        fs.append(optimal / ufloat(data.value, data.rel_unc if data.rel_unc != 0 else data.unc))
                        ds.append(optimal.s / data.unc)
                except KeyError:  # opt_rad_info doesn't have right form
                    pass
            if not fs:
                continue
            print(z, fs, ds)
            f_factors[z] = sum(f / f.s ** 2 for f in fs) / sum(1 / f.s**2 for f in fs)  # weighted average
            delta_factors[z] = ufloat(np.mean(ds), np.std(ds))
        return f_factors, delta_factors


class DataGrouper:

    @staticmethod
    def group(rad_info: RadiiInformation) -> list[DataGroup]:

        def build_graph() -> tuple[dict[Nuclide, set[Nuclide]], dict[Nuclide, set[RadiusData]]]:
            graph_edges = defaultdict(set)
            graph_nodes = defaultdict(set)

            for data in rad_info.values():
                nuclides = data.nuclides
                for nuc1, nuc2 in product(nuclides, nuclides):
                    if nuc1 == nuc2:
                        graph_nodes[nuc1].add(data)
                    else:
                        graph_edges[nuc1].add(nuc2)
                        graph_edges[nuc2].add(nuc1)
            return graph_edges, graph_nodes

        def graph_bfs() -> list[RadiiInformation]:
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

        def identify_primary_and_secondary(info: RadiiInformation) -> DataGroup:
            """Categorize data as being primary (will affect chi-squared) or secondary (will not affect chi-squared).
            First searches through nuclides to see which are primary/secondary, then assigns data as primary iff
            all nuclides are primary nuclides."""

            # absolute measurements make a nuclide primary
            primary_nuclides = set(itertools.chain.from_iterable([data.nuclides for data in info.values()
                                                                  if data.term.termtype == 'absolute']))

            if not primary_nuclides:
                return DataGroup(all_data=info, primary_data=info)

            secondary_nuclides = []
            # strategy is to simultaenous floodfill from all primary nuclides, each worker deleting edges behind them
            # if two workers meet, then all nuclides in the path are considered primary
            # if a worker has nowhere to go on its turn, and is on a primary nuclide, all nuclides in its path are primary
            # if a worker has nowhere to go on its turn, but is not on a primary nuclide, then any nuclides in its path that
            # are not primary nuclides MUST BE secondary nuclides, and we can add a list of secondary nuclides
            queue = deque([(nuc, {nuc}) for nuc in primary_nuclides])
            num_present = defaultdict(int)
            for nuc in primary_nuclides:
                num_present[nuc] += 1

            while queue:
                loc, seen = queue.popleft()

                had_edges = edges[loc] != set()
                for next_nuc in list(edges[loc]):
                    next_seen = seen.copy() | {next_nuc}
                    queue.append((next_nuc, next_seen))
                    num_present[next_nuc] += 1
                    # burn the bridge
                    edges[loc].remove(next_nuc)
                    if loc in edges[next_nuc]:
                        edges[next_nuc].remove(loc)

                if num_present[loc] > 1:  # we are not alone ...
                    primary_nuclides.update(seen)
                    continue

                num_present[loc] -= 1

                if not had_edges:  # nowhere to go, no kids sent out
                    if loc in primary_nuclides:  # got somewhere already hit by other workers? might be redundant
                        primary_nuclides.update()
                    else:
                        secondary_nuclides.append(seen)

            secondary_nuclides = [[nuc for nuc in secondary if nuc not in primary_nuclides] for secondary in
                                  secondary_nuclides]
            secondary_nuclides = [set(group) for group in secondary_nuclides if group]

            # If there is a branch in the secondary data, some secondary nuclides may be in multiple entries of
            # secondary_nuclides -> we need to merge entries that share nuclides!
            changed = True
            while changed:
                changed = False
                for i, seci in list(enumerate(secondary_nuclides)):
                    for j, secj in list(enumerate(secondary_nuclides)):
                        if j <= i:
                            continue
                        if seci.intersection(secj):
                            secondary_nuclides.append(secondary_nuclides.pop(j).union(secondary_nuclides.pop(i)))
                            changed = True
                            break
                    if changed:
                        break

            primary_data, secondary_data = RadiiInformation(), [RadiiInformation() for _ in range(len(secondary_nuclides))]

            for data in info.values():
                if all(nuc in primary_nuclides for nuc in data.nuclides):
                    primary_data[data.data_id] = data
                    continue
                for i, secondary_group in enumerate(secondary_nuclides):
                    if any(nuc in secondary_group for nuc in data.nuclides):
                        secondary_data[i][data.data_id] = data
                        break
                else:
                    raise ValueError(f'Data: {data} not in any group!')

            if any(len(secondary) == 0 for secondary in secondary_data):
                raise ValueError('Some secondary information does not have any data')

            return DataGroup(all_data=info, primary_data=primary_data, secondary_data=secondary_data)

        edges, nodes = build_graph()
        graph_groups = graph_bfs()
        return [identify_primary_and_secondary(rad_info) for rad_info in graph_groups]

    def group_by_lab(self, rad_info: RadiiInformation) -> list[RadiiInformation]:
        pass

    @staticmethod
    def group_by_element(rad_info: RadiiInformation) -> list[DataGroup]:
        z_groups = defaultdict(list)
        for data in rad_info.values():
            atomic_numbers = [nuc.z for nuc in data.nuclides]
            if len(set(atomic_numbers)) > 1:  # this measurement extends between elements
                continue
            z_groups[atomic_numbers[0]].append(data)

        groups = []
        for atomic_number, rad_data in z_groups.items():
            group = DataGroup(all_data=rad_info.split(rad_data), primary_data=rad_info.split(rad_data),
                              name=f'Z={atomic_number:03}')
            groups.append(group)
        groups = sorted(groups, key=lambda g: g.name)

        return groups

    @staticmethod
    def group_by_term(rad_info: RadiiInformation) -> list[DataGroup]:
        terms = set([data.term.termstr for data in rad_info.values()])

        rad_infos = {term: rad_info.split([data for data in rad_info.values() if data.term.termstr == term]) for term in
                     terms}

        return [DataGroup(all_data=info, primary_data=info, name=f'Term: {term}') for
                term, info in rad_infos.items()]


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
