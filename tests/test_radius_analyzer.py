from radbase.radius_analyzer import Term, Nuclide, RadiusData, RadiiInformation, DataGrouper, RadiusAnalyzer
from lmfit import Parameters
import numpy as np

params = Parameters()
params.add('R001001', 1)
params.add('R001002', 2)


def test_terms():
    absterm = Term('R001001')
    assert absterm.eval(params) == 1
    assert absterm.termtype == 'absolute'

    linrelterm = Term('R001002-R001001')
    assert linrelterm.eval(params) == 1
    assert linrelterm.termtype == 'linear relative'

    sqrelterm = Term('R001002**2-R001001**2')
    assert sqrelterm.eval(params) == 3
    assert sqrelterm.termtype == 'squared relative'
    assert sqrelterm.get_nuclides() == [Nuclide('R001002'), Nuclide('R001001')]


def test_radius_data():
    data = RadiusData('R001001', 0.9, 0.1, 0)
    assert np.isclose(data.residual(params), 1.0, 1e-4)
    assert data.nuclides == [Nuclide('R001001')]


def test_radii_information():
    rad_info = RadiiInformation()
    rad_info.add('R001001', 0, 1)
    rad_info.add('R001002-R001001', 0, 1)
    rad_info.add('R001003-R001001', 0, 1)
    rad_info.add('R001004**2-R001005**2', 0, 1)

    assert sorted(rad_info.nuclides, key=lambda n: n.z) == [Nuclide(f'R00100{i}') for i in range(1, 6)]


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


def test_radius_analyzer():
    pass
