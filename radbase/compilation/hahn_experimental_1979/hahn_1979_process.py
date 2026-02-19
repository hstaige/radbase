import csv

from uncertainties import ufloat_fromstr

from radbase.data_entry import (DataEntryInterface,
                                muonic_transition_energy_template)

values: dict = {
    'Reference': 'hahn_experimental_1979'
}

interface = DataEntryInterface(start_interface=False)

with open('./energies.txt', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    elements = reader[0][1:]
    for row in reader[1:]:
        transition, energies = row[0].strip(), row[1:]

        transition = {side: level for side, level in zip(['Upper', 'Lower'], transition.split('-'))}
        values['Transition'] = transition

        for nuclide, energy in zip(elements, energies):
            nuclide, energy = nuclide.strip(), energy.strip()
            values['Nuclide'] = nuclide
            first_par, last_par = energy.find('('), energy.find(')')
            energy = ufloat_fromstr(energy[:first_par] + energy[last_par + 1:])
            values['Energy [keV]'] = {'Value': energy.n,
                                      'Uncertainty': energy.s}

            interface.save_data(muonic_transition_energy_template, values, replacement_strategy='AlwaysReplace')
