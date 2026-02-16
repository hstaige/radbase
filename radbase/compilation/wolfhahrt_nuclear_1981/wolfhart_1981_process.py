import csv

from radbase.data_entry import (DataEntryInterface,
                                muonic_barret_theory_template,
                                muonic_transition_energy_difference_template,
                                muonic_transition_energy_template)

ediff_values: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Transition': {'Upper': '2p3/2', 'Lower': '1s1/2'},
    'Notes': 'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting. '
}

bvalues: dict = {
    'Reference': 'wolfhahrt_nuclear_1981',
    'Nuclear Polarization Method': 'Calculated',
    'Calculated Nuclear Polarization': [],
    'Notes': 'No theoretical uncertainties. #Exp_unc'
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
