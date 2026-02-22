import re

import pandas as pd
from scipy.constants import c, electron_volt, hbar, pi
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin

from radbase.data_entry import (
    DataEntryInterface, electron_scattering_cross_section_ratio_template)


def angle_acceptance(quantity_str):
    angle_acceptances = {'ᵃ': 0.9, 'ᵇ': 0.3, 'ᶜ': 0.36, 'ᵈ': 0.2}
    if quantity_str[-1] in angle_acceptances:
        return quantity_str[:-1], angle_acceptances[quantity_str[-1]]
    else:
        return quantity_str, 1.8


def convert_value(quantity_str):
    quantity_str = quantity_str.replace('−', '-')
    if len(matches := re.findall(r'\([^)^(]*\)', quantity_str)) == 2:
        value, exponent = matches
        value = ufloat_fromstr(value[1:-1]) * 10 ** int(exponent[1:-1])
    else:
        value = ufloat_fromstr(quantity_str)
    return value


df198MeV = pd.read_csv('./cross_section_ratios_198MeV.csv')
df299p5MeV = pd.read_csv('./cross_section_ratios_299p5MeV.csv')

interface = DataEntryInterface(start_interface=False)
csr_procs = electron_scattering_cross_section_ratio_template.proc_dict

for (df, nrg) in [(df198MeV, 198), (df299p5MeV, 299.5)]:
    for col in df.columns:
        if col == 'theta':
            continue

        for i, row in df.iterrows():

            if row[col] is None or row[col] in ['', '...']:
                continue

            csr_values = {}
            csr_values |= csr_procs['Reference'].process_data('heisenberg_electron_1972')

            theta = row['theta']
            csr_values |= csr_procs['theta [deg]'].process_data(theta)

            nrg_j = nrg * electron_volt * 10 ** 6
            p_i = nrg_j / c  # Ultra-relativistic
            q = 2 * p_i / hbar * sin(float(theta) / 2 * pi / 180) * 10 ** -15  # in units of fm^-1
            csr_values |= csr_procs['q [1/fm]'].process_data(str(q))

            csr_values |= csr_procs['Energy [MeV]'].process_data(str(nrg))

            csr_values |= csr_procs['Nuclide A'].process_data(col[2:6])
            csr_values |= csr_procs['Nuclide B'].process_data(col[7:11])

            dratio = row[col]
            dratio = convert_value(dratio) / 100  # convert from percent
            ratio = (1 - dratio) / (1 + dratio)

            csr_values |= csr_procs['Cross section ratio (A/B) [-]'].process_data(str(ratio))

            interface.save_data(electron_scattering_cross_section_ratio_template, csr_values,
                                replacement_strategy='AlwaysReplace')
