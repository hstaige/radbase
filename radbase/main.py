from collections import defaultdict

from radbase.database import AngeliDataMapper
from radbase.radius_analyzer import (DataGrouper, RadiiInformation,
                                     RadiusAnalyzer, compare_radii_information)
from radbase.startup.load_results import (angeli_abs_info, angeli_rem_info,
                                          angeli_remko_info,
                                          angeli_remnoniso_info,
                                          angeli_rfinal_info)

if __name__ == '__main__':

    def prev_params(group, prev_opt):
        params = analyzer.guess_params(group)
        for par in params:
            match = [data for data in prev_opt.values() if data.term.termstr == par]
            if match:
                params[par].value = match[0].value
        return params

    measurements_loc = r'../inputs/measurements.csv'
    excel_loc = r'../inputs/NuclearRadius_all_sheets_renamed.xlsx'

    loader = AngeliDataMapper(measurements_loc)
    analyzer = RadiusAnalyzer()
    grouper = DataGrouper()

    # print('Absolute')
    # rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'isoe', 'isom', 'noniso'])
    # rad_info = analyzer.adjust_uncertainties(rad_info)
    # opt_rad_info = RadiiInformation().join([analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    # compare_radii_information('../outputs/absl_comparison.xlsx', [angeli_abs_info(excel_loc), opt_rad_info],
    #                           prefixes=['angeli_abs', 'hunter_abs'])

    print('Noniso: R+')
    rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'isoe', 'isom'])
    rad_info = analyzer.adjust_uncertainties(rad_info)
    opt_rad_info = RadiiInformation().join(
        [analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    compare_radii_information('../outputs/noniso_comparison.xlsx', [angeli_remnoniso_info(excel_loc), opt_rad_info],
                              prefixes=['angeli_noniso', 'hunter_noniso'])
    non_iso_info = opt_rad_info
    for data in non_iso_info.values():  # delete correlations so values are are independent.
        data.correlations = defaultdict(dict)
    print(non_iso_info)

    # print('REM: R++')
    # rad_info = loader.load_data(bad_types=['ois', 'kalpha'])
    # rad_info = analyzer.adjust_uncertainties(rad_info)
    # opt_rad_info = RadiiInformation().join(
    #     [analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    #
    # compare_radii_information('../outputs/rem_comparison.xlsx', [angeli_rem_info(excel_loc), opt_rad_info],
    #                           prefixes=['angeli_rem', 'hunter_rem'])
    # full_rem = opt_rad_info

    # print('REM: no noniso')
    # rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'noniso', 'absl'])
    # rad_info = analyzer.adjust_uncertainties(rad_info)
    # rad_info = rad_info.join([non_iso_info])
    # groups = grouper.group(rad_info)
    # for group in groups:
    #     group = sorted(list(group.values()), key=lambda d: d.data_id)
    #     for data in group:
    #         print(data)
    #     print('\n')
    # print('Group lengths:', sorted([len(group) for group in groups], reverse=True))
    # opt_rad_info = RadiiInformation().join(
    #     [analyzer.optimize_radii(rad_group, params=prev_params(rad_group, opt_rad_info)) for rad_group in groups])
    #
    # compare_radii_information('../outputs/rem_without_nonisocomparison.xlsx', [angeli_rem_info(excel_loc), opt_rad_info],
    #                           prefixes=['angeli_rem', 'hunter_rem'])

    # compare_radii_information('../outputs/rem_compare_separated.xlsx', [full_rem, opt_rad_info],
    #                           prefixes=['REM', 'REM_separated'])

    print('Remko')
    rad_info = loader.load_data()
    rad_info = analyzer.adjust_uncertainties(rad_info)
    groups = grouper.group(rad_info)
    groups = sorted(groups, key=lambda g: len(g.all_data))
    print([len(group.all_data) for group in groups])
    for group in groups:
        print('Total Nuclides:', len(group.all_data.nuclides))
        print('Primary Nuclides:', len(group.primary_data.nuclides))
        print('Primary Measurements:')
        primary = sorted(list(group.primary_data.values()), key=lambda d: d.data_id)
        for data in primary:
            print(data)
        print('\n')
        print('Secondary Nuclides:', len(group.nuclides) - len(group.primary_data.nuclides))
        print('Secondary Measurements:')
        for i, secondary in enumerate(group.secondary_data):
            print(f'Secondary Group {i}')
            secondary = sorted(list(secondary.values()), key=lambda d: d.data_id)
            for data in secondary:
                print(data)
        print('\n')
    opt_rad_info = RadiiInformation().join(
        [analyzer.optimize_radii(rad_group, params=prev_params(rad_group, opt_rad_info)) for rad_group in groups])
    print('Optimization finished!')

    compare_radii_information('../outputs/remko.xlsx',
                              [angeli_remko_info(excel_loc), opt_rad_info],
                              prefixes=['angeli_remko', 'hunter_remko'])

    f_factors, delta_factors = {}, {}
    for group in groups:
        f_facts, d_facts = analyzer.calculate_f_factors(group.all_data, opt_rad_info)
        f_factors.update(f_facts), delta_factors.update(d_facts)

        not_ois_nuclides = set()
        for data in group.all_data.values():
            if getattr(data, 'measurement_type', 'ois') != 'ois':
                not_ois_nuclides.update(set(data.nuclides))

        for secondary in group.secondary_data:
            for data in secondary.values():
                if data.measurement_type == 'ois' and not all(nuc in not_ois_nuclides for nuc in data.nuclides):
                    z = data.nuclides[0].z
                    if z in f_factors:
                        data.unc = ((data.unc * delta_factors[z].n) ** 2 + (data.value * f_factors[z].s) ** 2) ** (1/2)
                        data.value *= f_factors[z].n
    print(f_factors)

    opt_rad_info = RadiiInformation().join(
        [analyzer.optimize_radii(rad_group, params=prev_params(rad_group, opt_rad_info)) for rad_group in groups])
    print('Optimization finished!')

    compare_radii_information('../outputs/rfinal.xlsx',
                              [angeli_rfinal_info(excel_loc), opt_rad_info],
                              prefixes=['angeli_remko_final', 'hunter_remko_final'])
