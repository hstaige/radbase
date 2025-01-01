from radbase.radius_analyzer import create_term, Nuclide, RadiusData, RadiiInformation, DataGrouper, RadiusAnalyzer, \
    AbsoluteTerm, LinearRelativeTerm, SquaredRelativeTerm
from lmfit import Parameters
from uncertainties import ufloat, UFloat, correlated_values_norm
import numpy as np


def test_terms():
    params = Parameters()
    params.add('R001001', 1)
    params.add('R001002', 2)

    uvars = params.create_uvars(covar=np.array([[0.1**2, 0], [0, 0.1**2]]))

    absterm = create_term('R001001')
    assert absterm.asteval(params) == 1
    assert absterm.eval(params) == 1
    assert absterm.calc_uvar(uvars).nominal_value == 1
    assert np.isclose(absterm.calc_uvar(uvars).std_dev, 0.1)
    assert absterm.termtype == 'absolute'
    assert isinstance(absterm, AbsoluteTerm)

    linrelterm = create_term('R001002-R001001')
    assert linrelterm.asteval(params) == 1
    assert linrelterm.eval(params) == 1
    assert linrelterm.calc_uvar(uvars).nominal_value == 1
    assert np.isclose(linrelterm.calc_uvar(uvars).std_dev, np.sqrt(0.1 ** 2 + 0.1 ** 2))
    assert linrelterm.termtype == 'linear_relative'
    assert isinstance(linrelterm, LinearRelativeTerm)

    sqrelterm = create_term('R001002**2-R001001**2')
    assert sqrelterm.asteval(params) == 3
    assert sqrelterm.eval(params) == 3
    assert sqrelterm.calc_uvar(uvars).nominal_value == 3
    assert np.isclose(sqrelterm.calc_uvar(uvars).std_dev, np.sqrt(0.2 ** 2 + 0.4 ** 2))
    assert sqrelterm.termtype == 'squared_relative'
    assert sqrelterm.get_nuclides() == [Nuclide('R001002'), Nuclide('R001001')]
    assert isinstance(sqrelterm, SquaredRelativeTerm)


def test_radius_data():
    params = Parameters()
    params.add('R001001', 1)
    params.add('R001002', 2)

    data = RadiusData('R001001', 0.9, 0.1, 0)
    assert np.isclose(data.residual(params), 1.0, 1e-4)
    assert data.nuclides == [Nuclide('R001001')]


def test_radii_information():
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 2)
    rad_info.add('R001003-R001001', 0, 3)
    rad_info.add('R001004**2-R001005**2', 0, 4)

    assert sorted(rad_info.nuclides, key=lambda n: n.z) == [Nuclide(f'R00100{i}') for i in range(1, 6)]

    rad_info.correlate(0, 1, 0.3)
    rad_info.correlate(1, 2, 0.4)
    rad_info.correlate(3,2, 0.5)

    corr_matrix = np.array([[1, 0.3, 0, 0], [0.3, 1, 0.4, 0], [0, 0.4, 1, 0.5], [0, 0, 0.5, 1]])
    assert np.allclose(rad_info.get_correlation_matrix(), corr_matrix)

    uncs = np.diag([1, 2, 3, 4])
    cov_matrix = uncs @ corr_matrix @ uncs
    assert np.allclose(rad_info.get_covariance_matrix(), cov_matrix)

    uvars = rad_info.create_uvars()
    assert len(uvars) == len(rad_info)
    assert all(isinstance(uvar, UFloat) for uvar in uvars)

    uvars = correlated_values_norm([(0, 1), (0, 2), (0, 3), (0, 4)], corr_matrix)
    terms = [data.term for data in rad_info._get_data_sorted()]
    rad_info_from_uvar = RadiiInformation.from_terms_and_uvars(terms, uvars)
    assert [np.isclose(rad_info_from_uvar[i].value, rad_info[i].value) for i in rad_info]
    assert [np.isclose(rad_info_from_uvar[i].unc, rad_info[i].unc) for i in rad_info]
    assert [rad_info_from_uvar[i].correlations == rad_info[i].correlations for i in rad_info]


def test_data_grouper():
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001003-R001001', 0, 1)
    rad_info.add('R001004**2-R001005**2', 0, 1)
    rad_info.add('R001004**2-R001006**2', 0, 1)
    rad_info.add('R002004**2-R002006**2', 0, 1)

    groups = DataGrouper.group(rad_info)
    assert all(isinstance(g, RadiiInformation) for g in groups)
    groups = sorted(groups, key=lambda x: len(x.keys()))
    assert [len(g.keys()) for g in groups] == [1, 2, 3]

    groups = DataGrouper.group_by_isotope(rad_info)
    assert all(isinstance(g, RadiiInformation) for g in groups)
    groups = sorted(groups, key=lambda x: len(x.keys()))
    assert [len(g.keys()) for g in groups] == [1, 5]