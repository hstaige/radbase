import csv
import json

from uncertainties import ufloat, ufloat_fromstr

from radbase.data_entry import (
    DataEntryInterface, charge_distribution_template,
    electron_scattering_cross_section_ratio_template,
    electron_scattering_cross_section_template, muonic_barrett_shift_template,
    muonic_barrett_theory_template,
    muonic_nuclear_polarization_calculation_template,
    muonic_qed_calculation_template,
    muonic_transition_energy_difference_template,
    muonic_transition_energy_template, radius_difference_template)

cd_procs = charge_distribution_template.proc_dict
np_procs = muonic_nuclear_polarization_calculation_template.proc_dict
qed_procs = muonic_qed_calculation_template.proc_dict
bar_procs = muonic_barrett_theory_template.proc_dict
bardiff_procs = muonic_barrett_shift_template.proc_dict
e_procs = muonic_transition_energy_template.proc_dict
ediff_procs = muonic_transition_energy_difference_template.proc_dict
raddiff_procs = radius_difference_template.proc_dict

reference = 'wohlfahrt_nuclear_1981'
transition = {'Upper': '2p3/2', 'Lower': '1s1/2'}
transition_or_level = [transition, '']

ediff_values = {}
ediff_values |= ediff_procs['Reference'].process_data(reference)
ediff_values |= ediff_procs['Transition\nor Level'].process_data(transition_or_level)
ediff_values |= ediff_procs['Notes'].process_data(
    'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting.')

e_values = {}
e_values |= e_procs['Reference'].process_data(reference)
e_values |= e_procs['Transition'].process_data(transition)
e_values |= e_procs['Notes'].process_data(
    'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting.')

np_values = {}
np_values |= np_procs['Reference'].process_data(reference)
np_values |= np_procs['Transition\nor Level'].process_data(transition_or_level)

qed_values = {}
qed_values |= qed_procs['Reference'].process_data(reference)
qed_values |= qed_procs['Transition\nor Level'].process_data(transition_or_level)

bar_values = {}
bar_values |= bar_procs['Reference'].process_data(reference)

bardiff_values = {}
bardiff_values |= bardiff_procs['Reference'].process_data(reference)

cd_values = {}
cd_values |= cd_procs['Reference'].process_data(reference)
cd_values |= cd_procs['Notes'].process_data('a fixed, c varied to reproduce 2p3/2-1s1/2 energies.')

raddiff_values = {}
raddiff_values |= raddiff_procs['Reference'].process_data(reference)

interface = DataEntryInterface(start_interface=False)

with open('./energy_diffs.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    ediff_keys = {}
    ediff_procs = muonic_transition_energy_difference_template.proc_dict
    for row in reader:
        nuclideA, nuclideB, energy_diff, unc = row

        nuclideA = nuclideA[2:] + nuclideA[:2]
        nuclideB = nuclideB[2:] + nuclideB[:2]
        ediff_values |= ediff_procs['Nuclide A'].process_data(nuclideA)
        ediff_values |= ediff_procs['Nuclide B'].process_data(nuclideB)

        ediff_values |= ediff_procs['Energy Difference (A-B) [keV]'].process_data(
            repr(ufloat(float(energy_diff), float(unc))))

        interface.save_data(muonic_transition_energy_difference_template, ediff_values,
                            replacement_strategy='AlwaysReplace')

        ediff_keys[(nuclideA, nuclideB)] = muonic_transition_energy_difference_template.data_key(ediff_values)

with open('barrett_info.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    rely_dict = {}

    for row in reader:
        nuclide, energy, np, qed, c, r, k, alpha, cz, rka, *_ = row

        nuclide = nuclide[2:] + nuclide[:2]
        energy = ufloat_fromstr(energy)
        e_values |= e_procs['Nuclide'].process_data(nuclide)
        e_values |= e_procs['Energy [keV]'].process_data(repr(energy))
        e_values |= e_procs['Notes'].process_data(
            'Determined from least squares fitting. HF splitting and isotopic impurity corrected. 2p3/2-1s from CoM of 2p-1s and theoretical FS splitting.')

        interface.save_data(muonic_transition_energy_template,
                            e_values, replacement_strategy='AlwaysReplace')

        np_values |= np_procs['Nuclide'].process_data(nuclide)
        np_values |= np_procs['Energy [keV]'].process_data(repr(ufloat(float(np), float(np) * 0.3)))
        np_values |= np_procs['Notes'].process_data('30% uncertainty assigned.')

        interface.save_data(muonic_nuclear_polarization_calculation_template,
                            np_values, replacement_strategy='AlwaysReplace')

        qed_values |= qed_procs['Nuclide'].process_data(nuclide)
        qed_values |= qed_procs['Energy [keV]'].process_data(repr(ufloat(float(qed), 0.02)))
        qed_values['Notes'] = 'Page 539, constant uncertainty'

        interface.save_data(muonic_qed_calculation_template,
                            qed_values, replacement_strategy='AlwaysReplace')

        bar_values |= bar_procs['Nuclide'].process_data(nuclide)
        cd_values |= cd_procs['Nuclide'].process_data(nuclide)

        e_key = muonic_transition_energy_template.data_key(e_values)

        muonic_np_key = muonic_nuclear_polarization_calculation_template.data_key(np_values)
        muonic_qed_key = muonic_qed_calculation_template.data_key(qed_values)
        rely_on = {key: True for key in [e_key, muonic_np_key, muonic_qed_key]}
        rely_dict[nuclide] = rely_on

        rka = ufloat_fromstr(rka)
        rka += ufloat(0, 0.003)

        bar_values |= bar_procs['Relies On'].process_data(rely_on)
        cd_values |= cd_procs['Relies On'].process_data(rely_on)

        bar_values |= bar_procs['Rka [fm]'].process_data(repr(rka))
        bar_values |= bar_procs['k [-]'].process_data(k)
        bar_values |= bar_procs['alpha [1/fm]'].process_data(alpha)
        bar_values |= bar_procs['Cz [fm/keV]'].process_data(cz)
        bar_values |= bar_procs['Notes'].process_data(
            'Statistical uncertainty and 0.003 fm uncertainty (page 538, Table IV caption) combined in quadrature.')

        cd_values |= cd_procs['Relies On'].process_data({key: True for key in [e_key, muonic_np_key, muonic_qed_key]})

        interface.save_data(muonic_barrett_theory_template,
                            bar_values, replacement_strategy='AlwaysReplace')

        cd_values |= cd_procs['Charge Distribution Parameters'].process_data(
            {'cd': 'TwoParameterFermi', 'cd_data': [c, '0.55']})
        cd_values |= cd_procs['Notes'].process_data('a fixed at 0.55 fm, see page 539')

        interface.save_data(charge_distribution_template,
                            cd_values, replacement_strategy='AlwaysReplace')

with open('moment_diffs.csv', 'r') as f:
    reader = list(csv.reader(f, delimiter=','))

    bardiff_values_orig = bardiff_values
    raddiff_values_orig = raddiff_values

    for row in reader:
        nuclide, rka_diff, rms_diff = row
        bardiff_values = bardiff_values_orig.copy()
        raddiff_values = raddiff_values_orig.copy()

        nuclideA, nuclideB = [nuc[2:] + nuc[:2] for nuc in nuclide.split('-')]

        bardiff_values |= bardiff_procs['Nuclide A'].process_data(nuclideA)
        bardiff_values |= bardiff_procs['Nuclide B'].process_data(nuclideB)

        raddiff_values |= raddiff_procs['Nuclide A'].process_data(nuclideA)
        raddiff_values |= raddiff_procs['Nuclide B'].process_data(nuclideB)

        bar_rely = {key: value for key, value in
                    ({ediff_keys[(nuclideA, nuclideB)]: True} | rely_dict[nuclideA] | rely_dict[nuclideB]).items() if
                    '_muonic_' not in key}  # exclude absolute transition energies from Barrett moment shift dependency.
        bardiff_values |= bardiff_procs['Relies On'].process_data(bar_rely)

        rka_diff = ufloat_fromstr(rka_diff) * 1e-3  # Table VI is in units of 10^-3
        rka_diff = rka_diff + ufloat(0, 0.001)  # See Table VI caption, NP uncertainty
        bardiff_values |= bardiff_procs['Rka (A-B) [fm]'].process_data(repr(rka_diff))
        bardiff_values |= bardiff_procs['k [-]'].process_data(None)
        bardiff_values |= bardiff_procs['alpha [1/fm]'].process_data(None)
        bardiff_values |= bardiff_procs['Cz [fm/keV]'].process_data(None)
        bardiff_values |= bardiff_procs['Notes'].process_data(
            'Includes the 0.001 fm uncertainty from NP (page 540, Table VI caption).')

        interface.save_data(muonic_barrett_shift_template,
                            bardiff_values, replacement_strategy='AlwaysReplace')
        barshift_key = muonic_barrett_shift_template.data_key(bardiff_values)


        def add_electron_scattering_ratio_links(nuc_a: str, nuc_b: str) -> tuple[dict[str, bool], str]:

            def get_keys(na, nb):
                with open(file_loc := rel_files[(na, nb)], 'r') as file:
                    data = list(json.load(file).keys())
                used_files.add(file_loc)
                return {key: True for key in data if na in key and nb in key}

            used_files = set()

            ratio_name = electron_scattering_cross_section_ratio_template.name
            frosch = f'../frosch_electron_1968/{ratio_name}.json'
            heisenberg = f'../heisenberg_electron_1972/{ratio_name}.json'
            lightbody = f'../lightbody_elastic_1983/{ratio_name}.json'
            rel_files = {('Ca48', 'Ca40'): frosch,
                         ('Ca44', 'Ca40'): frosch,
                         ('Ca42', 'Ca40'): frosch,
                         ('Ti48', 'Ca48'): frosch,
                         ('Ti48', 'Ti46'): heisenberg,
                         ('Ti50', 'Ti48'): heisenberg,
                         ('Cr52', 'Cr50'): lightbody,
                         ('Cr54', 'Cr52'): lightbody}
            rel_files = rel_files | {(nuc2, nuc1): ref for (nuc1, nuc2), ref in rel_files.items()}  # either order works

            if (nuc_a, nuc_b) in rel_files:  # one of the ratios listed in the refs
                keys = get_keys(nuc_a, nuc_b)
                return keys, ','.join(sorted(list(used_files)))

            # otherwise, ratios were chained to reach the desired ratio ex. Ca44/Ca42 = Ca44/Ca40 * Ca40/Ca42
            chains = {('Ca44', 'Ca42'): (('Ca44', 'Ca40'), ('Ca42', 'Ca40')),
                      ('Ti46', 'Ca44'): (('Ti48', 'Ti46'), ('Ti48', 'Ca48'), ('Ca48', 'Ca40'), ('Ca44', 'Ca40')),
                      ('Ti50', 'Ca48'): (('Ti50', 'Ti48'), ('Ti48', 'Ca48'))}

            if (nuc_a, nuc_b) in chains:
                result = {}
                for (nuc1, nuc2) in chains[(nuc_a, nuc_b)]:
                    result |= get_keys(nuc1, nuc2)
                return result, ','.join(sorted(list(used_files)))

            raise ValueError(f'{nuc_a}, {nuc_b} is not one of the pairs with ratios in Frosch68, Lightbody83, '
                             f'or Heisenberg1972. Are you sure this is a model independent combination?')


        raddiff_rely = {barshift_key: True}
        if rms_diff[-1] == 'a':
            rms_diff = rms_diff[:-1]
            ratio_keys, files_used = add_electron_scattering_ratio_links(nuclideA, nuclideB)
            raddiff_rely |= ratio_keys
            raddiff_values |= raddiff_procs['Relies On'].process_data(raddiff_rely)
            raddiff_values |= raddiff_procs['Notes'].process_data(
                f'Model independent radius difference determined from the cross section ratios in {files_used}.')
        else:
            raddiff_values |= raddiff_procs['Relies On'].process_data(raddiff_rely)
            raddiff_values |= raddiff_procs['Notes'].process_data('Model dependent, muonic 2pF only.')

        rms_diff = ufloat_fromstr(rms_diff) * 1e-3  # Table VI is in units of 10^-3
        rms_diff = rms_diff + ufloat(0, 0.001)  # See Table VI caption, NP uncertainty
        raddiff_values |= raddiff_procs['Radius Difference (A-B) [fm]'].process_data(repr(rms_diff))

        interface.save_data(radius_difference_template,
                            raddiff_values, replacement_strategy='AlwaysReplace')
