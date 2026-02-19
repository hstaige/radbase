import pandas as pd
from scipy.constants import c, electron_volt, hbar, pi
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin

from radbase.data_entry import (
    DataEntryInterface, electron_scattering_cross_section_ratio_template)

interface = DataEntryInterface(start_interface=False)

df = pd.read_csv('./cross_section_ratios.csv')

for i, row in df.iterrows():

    theta, nrg = row['theta'], row['Energy [MeV]']
    theta = ufloat(theta, 1)
    nrg_j = nrg * electron_volt * 10 ** 6
    p_i = nrg_j / c  # Ultra-relativistic
    q = 2 * p_i / hbar * sin(theta / 2 * pi / 180) * 10 ** -15  # in units of fm^-1

    for col in ['Cr52/Cr50', 'Cr54/Cr52', 'Cr54/Cr50']:
        cs_values: dict = {'Reference': 'lightbody_elastic_1983'}

        cs_values['q [1/fm]'] = {'Value': q.n,
                                 'Uncertainty': q.s}

        cs_values['Energy [MeV]'] = {'Value': nrg,
                                     'Uncertainty': None}
        cs_values['theta [deg]'] = {'Value': theta.n,
                                    'Uncertainty': theta.s}

        cs_values['Nuclide_A'] = col[:4]
        cs_values['Nuclide_B'] = col[5:]

        ratio = ufloat_fromstr(row[col])

        cs_values['Cross section Ratio (A/B) [fm^2/sr]'] = {'Value': ratio.n,
                                                            'Uncertainty': ratio.s}

        interface.save_data(electron_scattering_cross_section_ratio_template, cs_values,
                            replacement_strategy='AlwaysReplace')
