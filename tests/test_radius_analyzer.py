import numpy as np
from lmfit import Parameters
from uncertainties import UFloat, correlated_values_norm

from radbase.radius_analyzer import (AbsoluteTerm, DataGroup, DataGrouper,
                                     LinearRelativeTerm, Nuclide,
                                     RadiiInformation, RadiusAnalyzer,
                                     RadiusData, SquaredRelativeTerm,
                                     create_term)


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
    rad_info.correlate(3, 2, 0.5)

    corr_matrix = np.array([[1, 0.3, 0, 0], [0.3, 1, 0.4, 0], [0, 0.4, 1, 0.5], [0, 0, 0.5, 1]])
    assert np.allclose(rad_info.get_correlation_matrix(), corr_matrix)

    uncs = np.diag([1, 2, 3, 4])
    cov_matrix = uncs @ corr_matrix @ uncs
    assert np.allclose(rad_info.get_covariance_matrix(), cov_matrix)

    uvars = rad_info.create_uvars()
    assert len(uvars) == len(rad_info)
    assert all(isinstance(uvar, UFloat) for uvar in uvars)

    uvars = correlated_values_norm([(0, 1), (0, 2), (0, 3), (0, 4)], corr_matrix)
    terms = [data.term for data in rad_info.get_data_sorted()]
    rad_info_from_uvar = RadiiInformation.from_terms_and_uvars(terms, uvars)
    assert [np.isclose(rad_info_from_uvar[i].value, rad_info[i].value) for i in rad_info]
    assert [np.isclose(rad_info_from_uvar[i].unc, rad_info[i].unc) for i in rad_info]
    assert [rad_info_from_uvar[i].correlations == rad_info[i].correlations for i in rad_info]

    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002', 3, .1)
    rad_info.correlate(0, 1, 0.5)
    terms = [create_term('R001001'), create_term('R001002'), create_term('R001002-R001001')]
    rad_info_res = rad_info.evaluate_terms(terms)
    assert len(rad_info_res) == len(terms)
    assert np.allclose([data.value for data in rad_info_res.get_data_sorted()], [1, 3, 2])
    assert np.allclose([data.unc for data in rad_info_res.get_data_sorted()],
                       [0.1, 0.1, 0.1])

    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002', 3, .1)

    other_rad_info = RadiiInformation()
    other_rad_info.add('R001002-R001001', 2, 0.3)
    rad_info = RadiiInformation().join([rad_info, other_rad_info])
    assert len(rad_info) == 3
    assert [data.data_id == data_id for data_id, data in rad_info.items()]


def test_radius_analyzer():
    rad_analyzer = RadiusAnalyzer()

    # Test optimization with uncorrelated variables and len(rad_info.nuclides) == len(rad_info)
    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002-R001001', 3, .1)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert len(opt_rad_info) == len(rad_info.nuclides)
    assert np.allclose([opt_rad_info[0].value, opt_rad_info[1].value], [1, 4], rtol=1e-4)
    assert np.allclose([opt_rad_info[0].unc, opt_rad_info[1].unc], [0.1, 0.1 * np.sqrt(2)], rtol=1e-4)

    # Test optimization with uncorrelated variables and len(rad_info.nuclides) < len(rad_info)
    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002-R001001', 3, .1)
    rad_info.add('R001002', 4.3, .1)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert len(opt_rad_info) == len(rad_info.nuclides)
    assert np.allclose([opt_rad_info[0].value, opt_rad_info[1].value], [1.1, 4.2], rtol=1e-4)
    assert np.allclose([opt_rad_info[0].unc, opt_rad_info[1].unc], [0.08165, 0.08165], rtol=1e-4)

    # Test optimization with correlated variables and len(rad_info.nuclides) < len(rad_info)
    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002-R001001', 3, .1)
    rad_info.add('R001002', 4.3, .1)
    rad_info.correlate(0, 1, 0.5)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert len(opt_rad_info) == len(rad_info.nuclides)
    assert np.allclose([opt_rad_info[0].value, opt_rad_info[1].value], [1.1125, 4.225], rtol=1e-3)
    assert np.allclose([opt_rad_info[0].unc, opt_rad_info[1].unc], [0.06614, 0.0866], rtol=1e-3)
    diff_info = opt_rad_info.evaluate_terms(['R001002-R001001'])
    assert np.isclose(diff_info[0].value, 3.1125, rtol=1e-3)
    assert np.isclose(diff_info[0].unc, 0.06614, rtol=1e-3)

    # Test optimization with unccorrelated variables, linear relative term
    rad_info = RadiiInformation()
    rad_info.add('R001002-R001001', 3, .1)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert np.isclose(opt_rad_info[0].value, 3)
    assert np.isclose(opt_rad_info[0].unc, .1)

    # Test optimization with unccorrelated variables, squared relative term
    rad_info = RadiiInformation()
    rad_info.add('R001002**2-R001001**2', 3, .1)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert np.isclose(opt_rad_info[0].value, 3)
    assert np.isclose(opt_rad_info[0].unc, .1)

    # Test optimization with unccorrelated variables, squared relative term
    rad_info = RadiiInformation()
    rad_info.add('R001002-R001001', 1, .1)
    rad_info.add('R001002**2-R001001**2', 3, .1)
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert np.allclose([opt_rad_info[0].value, opt_rad_info[1].value], [1, 2], rtol=1e-3)
    assert np.allclose([opt_rad_info[0].unc, opt_rad_info[1].unc], [0.206, 0.1118], rtol=1e-3)

    # Test uncertainty adjustment by term
    rad_info = RadiiInformation()
    rad_info.add('R001001', 2.590, .0500)
    rad_info.add('R001001', 2.4920, .0111)

    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert np.isclose(opt_rad_info[0].unc, 0.010836, rtol=1e-3)
    adj_rad_info = rad_analyzer.adjust_uncertainties(rad_info)
    opt_rad_info = rad_analyzer.optimize_radii(adj_rad_info)
    assert np.isclose(opt_rad_info[0].unc, 0.01654, rtol=1e-3)

    # Test optimization with unccorrelated variables, squared relative term
    rad_info = RadiiInformation()
    rad_info.add('R001001', 1, .1)
    rad_info.add('R001002', 2, .1)
    rad_info.add('R001002**2-R001001**2', 2.5, 10)
    rad_info[0].measurement_type = 'mu'
    rad_info[1].measurement_type = 'mu;'
    rad_info[2].measurement_type = 'ois'
    opt_rad_info = rad_analyzer.optimize_radii(rad_info)
    assert np.isclose(opt_rad_info[0].value, 1, rtol=1e-3)
    assert np.isclose(opt_rad_info[0].unc, .1, rtol=1e-3)
    assert np.isclose(opt_rad_info[1].value, 2, rtol=1e-3)
    assert np.isclose(opt_rad_info[1].unc, .1, rtol=1e-3)
    assert np.isclose(rad_analyzer.calculate_f_factors(rad_info, opt_rad_info)[1].n, 3 / 2.5, rtol=1e-3)


def test_data_grouper():
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001003-R001001', 0, 1)
    rad_info.add('R001004**2-R001005**2', 0, 1)
    rad_info.add('R001004**2-R001006**2', 0, 1)
    rad_info.add('R002004**2-R002006**2', 0, 1)

    groups = DataGrouper.group(rad_info)
    assert all(isinstance(g, DataGroup) for g in groups)
    groups = sorted(groups, key=lambda x: len(x.all_data))
    assert [len(g.all_data) for g in groups] == [1, 2, 4]
    assert [len(g.primary_data) for g in groups] == [1, 2, 1]
    assert [len(g.secondary_data) for g in groups] == [0, 0, 2]

    groups = DataGrouper.group_by_element(rad_info)
    assert all(isinstance(g, DataGroup) for g in groups)
    groups = sorted(groups, key=lambda x: len(x.all_data))
    assert [len(g.all_data) for g in groups] == [1, 6]
    assert [len(g.primary_data) for g in groups] == [1, 6]

    groups = DataGrouper.group_by_term(rad_info)
    assert all(isinstance(g, DataGroup) for g in groups)
    groups = sorted(groups, key=lambda x: len(x.all_data))
    assert [len(g.all_data) for g in groups] == [1, 1, 1, 1, 1, 2]
    assert [len(g.primary_data) for g in groups] == [1, 1, 1, 1, 1, 2]
    assert [len(g.secondary_data) for g in groups] == [0, 0, 0, 0, 0, 0]


def test_grouping_edge_cases():
    # test no loop with one absolute
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001002-R001003', 0, 1)

    group = DataGrouper.group(rad_info)[0]
    assert len(group.all_data) == 3
    assert len(group.primary_data) == 1
    assert len(group.secondary_data) == 1
    assert len(group.secondary_data[0]) == 2

    # test loop with one absolute
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001002-R001003', 0, 1)
    rad_info.add('R001001-R001003', 0, 1)

    group = DataGrouper.group(rad_info)[0]
    assert len(group.all_data) == 4
    assert len(group.primary_data) == 4
    assert len(group.secondary_data) == 0

    # test branching with one absolute
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001002-R001003', 0, 1)
    rad_info.add('R001002-R001004', 0, 1)

    group = DataGrouper.group(rad_info)[0]
    assert len(group.all_data) == 4
    assert len(group.primary_data) == 1
    assert len(group.secondary_data) == 1
    assert len(group.secondary_data[0]) == 3

    # test branching with loop
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001002-R001003', 0, 1)
    rad_info.add('R001002-R001004', 0, 1)
    rad_info.add('R001004**2-R001003**2', 0, 1)

    group = DataGrouper.group(rad_info)[0]
    assert len(group.all_data) == 5
    assert len(group.primary_data) == 5
    assert len(group.secondary_data) == 0
