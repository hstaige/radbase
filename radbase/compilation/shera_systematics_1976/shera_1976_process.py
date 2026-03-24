import csv
import json

import numpy as np
from uncertainties import ufloat, ufloat_fromstr

from radbase.data_entry import (
    DataEntryInterface, charge_distribution_template,
    muonic_barrett_shift_template, muonic_barrett_theory_template,
    muonic_nuclear_polarization_calculation_template,
    muonic_qed_calculation_template,
    muonic_transition_energy_difference_template,
    muonic_transition_energy_template)

cd_procs = charge_distribution_template.proc_dict
np_procs = muonic_nuclear_polarization_calculation_template.proc_dict
qed_procs = muonic_qed_calculation_template.proc_dict
bar_procs = muonic_barrett_theory_template.proc_dict
bardiff_procs = muonic_barrett_shift_template.proc_dict
e_procs = muonic_transition_energy_template.proc_dict
ediff_procs = muonic_transition_energy_difference_template.proc_dict

reference = 'shera_systematics_1976'
transition_ka1 = {'Upper': '2p3/2', 'Lower': '1s1/2'}
transition_ka2 = {'Upper': '2p1/2', 'Lower': '1s1/2'}
transition_or_level_ka1 = [transition_ka1, '']
transition_or_level_ka2 = [transition_ka2, '']

ediff_values = {}
ediff_values |= ediff_procs['Reference'].process_data(reference)
ediff_values |= ediff_procs['Transition\nor Level'].process_data(transition_or_level_ka1)
ediff_values |= ediff_procs['Notes'].process_data('Determined from least squares fitting. HF splitting and isotopic impurity corrected.')

e_values = {}
e_values |= e_procs['Reference'].process_data(reference)
e_values |= e_procs['Transition'].process_data(transition_ka1)

np_values = {}
np_values |= np_procs['Reference'].process_data(reference)
np_values |= np_procs['Transition\nor Level'].process_data(transition_or_level_ka1)

qed_values = {}
qed_values |= qed_procs['Reference'].process_data(reference)
qed_values |= qed_procs['Transition\nor Level'].process_data(transition_or_level_ka1)

bar_values = {}
bar_values |= bar_procs['Reference'].process_data(reference)

bardiff_values = {}
bardiff_values |= bardiff_procs['Reference'].process_data(reference)

cd_values = {}
cd_values |= cd_procs['Reference'].process_data(reference)
cd_values |= cd_procs['Notes'].process_data('a fixed, c varied to reproduce 2p3/2-1s1/2 energies.')

interface = DataEntryInterface(start_interface=False)

combined_keys = {}
with open('energies.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    for row in reader:
        nuclide, ka2_nrg, ka1_nrg, combined_nrg = row

        nuclide = nuclide[2:] + nuclide[:2]
        e_values |= e_procs['Nuclide'].process_data(nuclide)

        e_values |= e_procs['Transition'].process_data(transition_ka2)
        e_values |= e_procs['Energy [keV]'].process_data(ka2_nrg)
        e_values |= e_procs['Notes'].process_data('')
        interface.save_data(muonic_transition_energy_template, e_values,
                            replacement_strategy='AlwaysReplace')

        e_values |= e_procs['Transition'].process_data(transition_ka1)
        e_values |= e_procs['Energy [keV]'].process_data(ka1_nrg)
        interface.save_data(muonic_transition_energy_template, e_values,
                            replacement_strategy='AlwaysReplace')

        e_values |= e_procs['Transition'].process_data(transition_ka1)
        e_values |= e_procs['Energy [keV]'].process_data(combined_nrg)
        e_values |= e_procs['Notes'].process_data('Combined Ka1 and Ka2 using theoretical FS splitting.')
        combined_keys[nuclide] = muonic_transition_energy_template.data_key(e_values)

        interface.save_data(muonic_transition_energy_template, e_values,
                            replacement_strategy='Suffix')


ediff_keys = {}
ediff_procs = muonic_transition_energy_difference_template.proc_dict
with open('./energy_diffs.csv', 'r') as f:
    data = f.read().splitlines()
    data = [line.split(',') for line in data]

    col_nucs, row_nucs = data[0][1:], [line[0] for line in data[1::2]]
    data = [[f'{r1}({r2})' for r1, r2 in zip(line1[1:], line2[1:])] for line1, line2 in zip(data[1::2], data[2::2])]  # merge uncertainties and values
    for i, row_nuc in enumerate(row_nucs):
        for j, col_nuc in enumerate(col_nucs):
            if j > i:
                continue
            e_diff = data[i][j]

            nuclideA = col_nuc[2:] + col_nuc[:2]
            nuclideB = row_nuc[2:] + row_nuc[:2]
            ediff_values |= ediff_procs['Nuclide A'].process_data(nuclideA)
            ediff_values |= ediff_procs['Nuclide B'].process_data(nuclideB)

            ediff_values |= ediff_procs['Energy Difference (A-B) [keV]'].process_data(e_diff)

            interface.save_data(muonic_transition_energy_difference_template, ediff_values,
                                replacement_strategy='AlwaysReplace')

            ediff_keys[(nuclideA, nuclideB)] = muonic_transition_energy_difference_template.data_key(ediff_values)
            ediff_keys[(nuclideB, nuclideA)] = muonic_transition_energy_difference_template.data_key(ediff_values)

with open('barrett.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    rely_dict = {}

    for row in reader:
        nuclide, energy, nucpol, qed, c, r, alpha, k, cz, rka, *_ = row

        cz = float(cz) / 1000

        nuclide = nuclide[2:] + nuclide[:2]

        interface.save_data(muonic_transition_energy_template,
                            e_values, replacement_strategy='AlwaysReplace')

        np_values |= np_procs['Nuclide'].process_data(nuclide)
        np_unc = float(nucpol) * 0.4
        np_values |= np_procs['Energy [keV]'].process_data(repr(ufloat(float(nucpol), np_unc)))
        np_values |= np_procs['Notes'].process_data('40% uncertainty assigned, see Page 736.')

        interface.save_data(muonic_nuclear_polarization_calculation_template,
                            np_values, replacement_strategy='AlwaysReplace')

        qed_values |= qed_procs['Nuclide'].process_data(nuclide)

        element_to_z = {'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30}
        zs = np.array([26, 28, 30])
        nuclear_recoils = np.array([0.055, 0.065, 0.077])
        qed_uncs = np.array([0.079, 0.091, 0.101])
        qed_val = float(qed) + np.interp(element_to_z[nuclide[:2]], zs, nuclear_recoils)
        qed_unc = np.interp(element_to_z[nuclide[:2]], zs, qed_uncs)

        qed_values |= qed_procs['Energy [keV]'].process_data(repr(ufloat(qed_val, qed_unc)))
        qed_values |= qed_procs['Notes'].process_data('Table VII \'all other corrections\' minus interpolated nuclear recoil. Uncertainty estimated from Table I.')

        interface.save_data(muonic_qed_calculation_template,
                            qed_values, replacement_strategy='AlwaysReplace')

        bar_values |= bar_procs['Nuclide'].process_data(nuclide)
        cd_values |= cd_procs['Nuclide'].process_data(nuclide)

        e_key = combined_keys[nuclide]
        muonic_np_key = muonic_nuclear_polarization_calculation_template.data_key(np_values)
        muonic_qed_key = muonic_qed_calculation_template.data_key(qed_values)
        rely_on = {key: True for key in [e_key, muonic_np_key, muonic_qed_key]}
        rely_dict[nuclide] = rely_on

        rka = ufloat_fromstr(rka)
        rka += ufloat(0, qed_unc * abs(float(cz)) / 1000)
        rka += ufloat(0, np_unc * abs(float(cz)) / 1000)

        bar_values |= bar_procs['Relies On'].process_data(rely_on)
        cd_values |= cd_procs['Relies On'].process_data(rely_on)

        bar_values |= bar_procs['Rka [fm]'].process_data(repr(rka))
        bar_values |= bar_procs['k [-]'].process_data(k)
        bar_values |= bar_procs['alpha [1/fm]'].process_data(alpha)
        bar_values |= bar_procs['Cz [fm/keV]'].process_data(float(cz) / 1000)
        bar_values |= bar_procs['Notes'].process_data(
            'Statistical uncertainty, NP uncertainty, and QED uncertainty combined in quadrature.')

        cd_values |= cd_procs['Relies On'].process_data({key: True for key in [e_key, muonic_np_key, muonic_qed_key]})

        interface.save_data(muonic_barrett_theory_template,
                            bar_values, replacement_strategy='AlwaysReplace')

        cd_values |= cd_procs['Charge Distribution Parameters'].process_data(
            {'cd': 'TwoParameterFermi', 'cd_data': [c, '0.55']})
        cd_values |= cd_procs['Notes'].process_data('a fixed at 0.55 fm, see page 743')

        interface.save_data(charge_distribution_template,
                            cd_values, replacement_strategy='AlwaysReplace')

with open('barrett_diffs.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    bardiff_values_orig = bardiff_values

    for row in reader:
        nuclide, rka_diff = row
        bardiff_values = bardiff_values_orig.copy()

        nuclideA, nuclideB = [nuc[2:] + nuc[:2] for nuc in nuclide.split('-')]

        bardiff_values |= bardiff_procs['Nuclide A'].process_data(nuclideA)
        bardiff_values |= bardiff_procs['Nuclide B'].process_data(nuclideB)

        bar_rely = {key: value for key, value in
                    ({ediff_keys[(nuclideA, nuclideB)]: True} | rely_dict[nuclideA] | rely_dict[nuclideB]).items() if
                    '_muonic_' not in key}  # exclude absolute transition energies from Barrett moment shift dependency.
        bardiff_values |= bardiff_procs['Relies On'].process_data(bar_rely)

        rka_diff = ufloat_fromstr(rka_diff) * 1e-3  # Table VI is in units of 10^-3
        bardiff_values |= bardiff_procs['Rka (A-B) [fm]'].process_data(repr(rka_diff))
        bardiff_values |= bardiff_procs['k [-]'].process_data(None)
        bardiff_values |= bardiff_procs['alpha [1/fm]'].process_data(None)
        bardiff_values |= bardiff_procs['Cz [fm/keV]'].process_data(None)
        bardiff_values |= bardiff_procs['Notes'].process_data(
            '#statonly')

        interface.save_data(muonic_barrett_shift_template,
                            bardiff_values, replacement_strategy='AlwaysReplace')
