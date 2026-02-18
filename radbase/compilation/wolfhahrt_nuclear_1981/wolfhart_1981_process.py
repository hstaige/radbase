import csv

from uncertainties import ufloat_fromstr

from radbase.data_entry import (
    DataEntryInterface, muonic_barret_theory_template,
    muonic_fermi_distribution_template,
    muonic_nuclear_polarization_calculation_template,
    muonic_qed_calculation_template,
    muonic_transition_energy_difference_template,
    muonic_transition_energy_template)

ediff_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Transition': {'Upper': '2p3/2', 'Lower': '1s1/2'},
    'Notes': 'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting. '
}

e_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Transition': {'Upper': '2p3/2', 'Lower': '1s1/2'},
    'Notes': 'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting. '
}

np_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Transition': {'Upper': '2p3/2', 'Lower': '1s1/2'},
}

qed_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Transition': {'Upper': '2p3/2', 'Lower': '1s1/2'},
}

b_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Nuclear Polarization Method': 'Calculated',
    'Notes': 'No NP or QED uncertainties included. #Exp_unc'
}

f_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Nuclear Polarization Method': 'Calculated',
    'a [fm]': {'Value': 0.55,
               'Uncertainty': None},
    'Notes': 'a fixed, c varied to reproduce 2p3/2-1s1/2 energies.'
}

interface = DataEntryInterface(start_interface=False)

with open('./energy_diffs.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    for row in reader:
        nuclideA, nuclideB, energy_diff, unc = row

        nuclideA = nuclideA[2:] + nuclideA[:2]
        nuclideB = nuclideB[2:] + nuclideB[:2]
        ediff_values['Nuclide_A'], ediff_values['Nuclide_B'] = nuclideA, nuclideB

        ediff_values['Energy Difference [keV] (A-B)'] = {'Value': energy_diff,
                                                         'Uncertainty': unc}

        interface.save_data(muonic_transition_energy_difference_template, ediff_values, replacement_strategy='AlwaysReplace')

with open('barrett_info.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    for row in reader:
        nuclide, energy, np, qed, c, r, k, alpha, cz, rka, *_ = row

        nuclide = nuclide[2:] + nuclide[:2]
        energy = ufloat_fromstr(energy)
        e_values['Nuclide'] = nuclide
        e_values['Energy [keV]'] = {'Value': energy.n,
                                    'Uncertainty': energy.s}

        interface.save_data(muonic_transition_energy_template,
                            e_values, replacement_strategy='AlwaysReplace')

        np_values['Nuclide'] = nuclide
        np_values['Energy [keV]'] = {'Value': float(np),
                                     'Uncertainty': float(np)*0.3}
        np_values['Notes'] = '30% uncertainty assigned.'

        interface.save_data(muonic_nuclear_polarization_calculation_template,
                            np_values, replacement_strategy='AlwaysReplace')

        qed_values['Nuclide'] = nuclide
        qed_values['Energy [keV]'] = {'Value': float(qed),
                                      'Uncertainty': float(qed)*0.3}
        qed_values['Notes'] = '30% uncertainty assigned.'

        interface.save_data(muonic_qed_calculation_template,
                            qed_values, replacement_strategy='AlwaysReplace')

        b_values['Nuclide'] = nuclide
        f_values['Nuclide'] = nuclide

        e_key = muonic_transition_energy_template.data_key(e_values)
        b_values['Previous Muonic Measurements'] = [e_key]
        f_values['Previous Muonic Measurements'] = [e_key]

        rka = ufloat_fromstr(rka)
        b_values['Rka [fm]'] = {'Value': rka.n,
                                'Uncertainty': rka.s}
        b_values['k [-]'] = k
        b_values['alpha [1/fm]'] = alpha
        b_values['Cz [fm/keV]'] = cz

        muonic_np_key = muonic_nuclear_polarization_calculation_template.data_key(np_values)
        b_values['Calculated Nuclear Polarization'] = [muonic_np_key]
        f_values['Calculated Nuclear Polarization'] = [muonic_np_key]

        muonic_qed_key = muonic_qed_calculation_template.data_key(qed_values)
        b_values['Calculated QED'] = [muonic_qed_key]
        f_values['Calculated QED'] = [muonic_qed_key]

        interface.save_data(muonic_barret_theory_template,
                            b_values, replacement_strategy='AlwaysReplace')

        f_values['c [fm]'] = {'Value': float(c),
                              'Uncertainty': None}

        interface.save_data(muonic_fermi_distribution_template,
                            f_values, replacement_strategy='AlwaysReplace')
