import re

import pandas as pd
from scipy.constants import c, electron_volt, hbar, pi
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin

from radbase.data_entry import (
    DataEntryInterface, electron_scattering_cross_section_ratio_template,
    electron_scattering_cross_section_template)

interface = DataEntryInterface(start_interface=False)
csr_procs = electron_scattering_cross_section_ratio_template.proc_dict
cs_procs = electron_scattering_cross_section_template.proc_dict

df = pd.read_csv('./cross_section_ratios.csv')

for i, row in df.iterrows():

    theta, nrg = row['theta'], row['Energy [MeV]']
    nrg_j = nrg * electron_volt * 10 ** 6
    p_i = nrg_j / c  # Ultra-relativistic
    q = 2 * p_i / hbar * sin(theta / 2 * pi / 180) * 10 ** -15  # in units of fm^-1

    for col in ['Cr52/Cr50', 'Cr54/Cr52', 'Cr54/Cr50']:
        csr_values = {}

        csr_values |= csr_procs['Reference'].process_data('lightbody_elastic_1983')
        csr_values |= csr_procs['q [1/fm]'].process_data(str(q))
        csr_values |= csr_procs['Energy [MeV]'].process_data(str(nrg))
        csr_values |= csr_procs['theta [deg]'].process_data(row['theta'])

        csr_values |= csr_procs['Nuclide A'].process_data(col[:4])
        csr_values |= csr_procs['Nuclide B'].process_data(col[5:])

        csr_values |= csr_procs['Cross section ratio (A/B) [-]'].process_data(row[col])

        interface.save_data(electron_scattering_cross_section_ratio_template, csr_values,
                            replacement_strategy='AlwaysReplace')


def convert_value(quantity_str):
    quantity_str = quantity_str.replace('−', '-')
    if len(matches := re.findall(r'\([^)^(]*\)', quantity_str)) == 2:
        value, exponent = matches
        value = ufloat_fromstr(value[1:-1]) * 10 ** int(exponent[1:-1])
    else:
        value = ufloat_fromstr(quantity_str[:-2]) * 10 ** int(quantity_str[-2:])
    return value


for nuclide, df in [('Cr50', pd.read_csv('./cross_section_Cr50.csv')),
                    ('Cr52', pd.read_csv('./cross_section_Cr52.csv')),
                    ('Cr54', pd.read_csv('./cross_section_Cr54.csv'))]:

    df.ffill(inplace=True)

    for i, row in df.iterrows():

        theta, nrg = row['theta'], row['Energy [MeV]']
        theta = ufloat(theta, 1)
        nrg_j = nrg * electron_volt * 10 ** 6
        p_i = nrg_j / c  # Ultra-relativistic
        q = 2 * p_i / hbar * sin(theta / 2 * pi / 180) * 10 ** -15  # in units of fm^-1

        cs_values = {}

        cs_values |= cs_procs['Reference'].process_data('lightbody_elastic_1983')
        cs_values |= cs_procs['q [1/fm]'].process_data(str(q))
        cs_values |= cs_procs['Energy [MeV]'].process_data(str(nrg))
        cs_values |= cs_procs['theta [deg]'].process_data(row['theta'])
        cs_values |= cs_procs['Nuclide'].process_data(nuclide)

        cs = convert_value(row['cross_section']) * 10  # mb to fm^2
        cs_values |= cs_procs['Cross section [fm^2/sr]'].process_data(str(cs))

        interface.save_data(electron_scattering_cross_section_template, cs_values,
                            replacement_strategy='AlwaysReplace')
