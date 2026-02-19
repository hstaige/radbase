import csv

from uncertainties import ufloat_fromstr

from radbase.data_entry import (DataEntryInterface,
                                muonic_barret_theory_template,
                                muonic_transition_energy_template)

evalues: dict = {
    'Reference': 'saito_muonic_2025'
}

bvalues: dict = {
    'Reference': 'saito_muonic_2025',
    'Nuclear Polarization Method': 'Calculated',
    'Calculated Nuclear Polarization': [],
    'Notes': 'No theoretical uncertainties. #Exp_unc'
}

interface = DataEntryInterface(start_interface=False)

with open('./values.txt', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    for row in reader:
        nuclide, transition, energy, k, alpha, rka = row

        nuclide = nuclide[3:] + nuclide[:-2]

        transition = {side: level for side, level in zip(['Upper', 'Lower'], transition.split('-'))}
        evalues['Transition'] = transition

        nuclide, energy = nuclide.strip(), energy.strip()
        evalues['Nuclide'], bvalues['Nuclide'] = nuclide, nuclide

        energy = ufloat_fromstr(energy)
        evalues['Energy [keV]'] = {'Value': energy.n,
                                   'Uncertainty': energy.s}

        interface.save_data(muonic_transition_energy_template, evalues, replacement_strategy='AlwaysReplace')

        bvalues['Previous Muonic Measurements'] = [muonic_transition_energy_template.data_key(evalues)]
        bvalues['alpha [1/fm]'] = alpha
        bvalues["k [-]"] = k
        rka = ufloat_fromstr(rka)
        bvalues['Rka [fm]'] = {'Value': rka.n,
                               'Uncertainty': rka.s}

        interface.save_data(muonic_barret_theory_template, bvalues, replacement_strategy='AlwaysReplace')
