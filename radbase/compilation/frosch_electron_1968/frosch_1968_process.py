import re
from pathlib import Path

import pandas as pd
from scipy.constants import c, electron_volt, hbar, pi
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin

from radbase.data_entry import (
    DataEntryInterface, electron_scattering_cross_section_ratio_template,
    electron_scattering_cross_section_template)


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


dfa = pd.read_csv('./cross_sections_a.csv')
dfb = pd.read_csv('./cross_sections_b.csv')
dfc = pd.read_csv('./cross_sections_c.csv')

for_all_cols = ['Energy [MeV]', 'theta']

by_quant_dfs = []
for df in [dfa, dfb, dfc]:
    for column in df.columns:
        if column not in for_all_cols:
            col_df = df[~(df[column].isin(['', '...',]) | df[column].isna())]
            col_df = col_df[for_all_cols + [column]]
            col_df.rename(columns={column: 'Value'}, inplace=True)
            col_df['Quantity'] = column
            by_quant_dfs.append(col_df)

by_quant_dfs = pd.concat(by_quant_dfs)
by_quant_dfs.to_csv('./cross_sections.csv', index=False)

# Remove previous data as we will be adding suffixes instead of replacing old data
Path(f'./{electron_scattering_cross_section_template.name}.json').unlink(missing_ok=True)
Path(f'./{electron_scattering_cross_section_ratio_template.name}.json').unlink(missing_ok=True)

interface = DataEntryInterface(start_interface=False)
cs_procs = electron_scattering_cross_section_template.proc_dict
csr_procs = electron_scattering_cross_section_ratio_template.proc_dict

for i, row in by_quant_dfs.iterrows():

    if row['Value'] is None:
        continue

    cs_values = {}
    cs_values |= cs_procs['Reference'].process_data('frosch_electron_1968')

    theta, nrg = row['theta'], row['Energy [MeV]']
    theta = ufloat(theta, 1)
    nrg_j = nrg * electron_volt * 10 ** 6
    p_i = nrg_j / c  # Ultra-relativistic
    q = 2 * p_i / hbar * sin(theta / 2 * pi / 180) * 10 ** -15  # in units of fm^-1

    cs_values |= cs_procs['q [1/fm]'].process_data(str(q))

    cs_values |= cs_procs['Energy [MeV]'].process_data(repr(nrg))
    cs_values |= cs_procs['theta [deg]'].process_data(repr(theta))

    if 'D' in row['Quantity']:
        cs_values |= csr_procs['Nuclide A'].process_data(row['Quantity'][2:6])
        cs_values |= csr_procs['Nuclide B'].process_data(row['Quantity'][7:11])
        cs_values |= csr_procs['q [1/fm]'].process_data(repr(q))

        dratio = row['Value']
        dratio, angular_acceptance = angle_acceptance(dratio)
        dratio = convert_value(dratio) / 100  # convert from percent
        ratio = (1 - dratio) / (1 + dratio)

        cs_values |= csr_procs['Cross section ratio (A/B) [-]'].process_data(repr(ratio))
        cs_values |= cs_procs['Notes'].process_data(f'Angular acceptance: {angular_acceptance:0.2f}')

        interface.save_data(electron_scattering_cross_section_ratio_template, cs_values,
                            replacement_strategy='Suffix')

    else:
        cs_values |= cs_procs['Nuclide'].process_data(row['Quantity'])

        cs = row['Value']
        cs, angular_acceptance = angle_acceptance(cs)
        cs = convert_value(cs)
        cs /= 10000  # microbarn to fm^2 conversion
        cs_values |= cs_procs['Cross section [fm^2/sr]'].process_data(repr(cs))

        cs_values |= cs_procs['Notes'].process_data(f'Angular acceptance: {angular_acceptance:0.2f}')

        interface.save_data(electron_scattering_cross_section_template, cs_values,
                            replacement_strategy='Suffix')
