import radbase
import pandas as pd


def test_groups_one_absolute():

    test_df = pd.DataFrame({'Iso_Idxs': ['001001']})
    answer = {'m_idxs': [0], 'iso_idxs': {'001001'}}

    group = radbase.solve.find_groups(test_df)[0]
    assert group['m_idxs'] == answer['m_idxs']
    assert group['iso_idxs'] == answer['iso_idxs']


def test_groups_one_relative():

    test_df = pd.DataFrame({'Iso_Idxs': ['001001,001002']})
    answer = {'m_idxs': [0], 'iso_idxs': {'001001', '001002'}}

    group = radbase.solve.find_groups(test_df)[0]
    assert group['m_idxs'] == answer['m_idxs']
    assert group['iso_idxs'] == answer['iso_idxs']


def test_groups_one_of_each():

    test_df = pd.DataFrame({'Iso_Idxs': ['001001,001002', '001001']})
    answer = {'m_idxs': [0, 1], 'iso_idxs': {'001001', '001002'}}

    group = radbase.solve.find_groups(test_df)[0]
    assert group['m_idxs'] == answer['m_idxs']
    assert group['iso_idxs'] == answer['iso_idxs']


def test_groups_seperate():

    test_df = pd.DataFrame({'Iso_Idxs': ['001001,001002', '001003']})
    answer1 = {'m_idxs': [1], 'iso_idxs': {'001003'}}
    answer2 = {'m_idxs': [0], 'iso_idxs': {'001001', '001002'}}

    groups = radbase.solve.find_groups(test_df)

    groups = sorted(groups, key=lambda g: len(g['iso_idxs']))

    group = groups[0]
    assert group['m_idxs'] == answer1['m_idxs']
    assert group['iso_idxs'] == answer1['iso_idxs']
    group = groups[1]
    assert group['m_idxs'] == answer2['m_idxs']
    assert group['iso_idxs'] == answer2['iso_idxs']