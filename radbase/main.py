from radbase.radius_analyzer import RadiiInformation, RadiusAnalyzer, DataGrouper, compare_radii_information
from radbase.database import AngeliDataMapper
from radbase.startup.load_results import angeli_abs_info, angeli_remnoniso_info, angeli_rem_info
from collections import defaultdict

if __name__ == '__main__':

    def prev_params(group, prev_opt):
        params = analyzer.guess_params(group)
        for par in params:
            match = [data for data in prev_opt.values() if data.term.termstr == par]
            if match:
                params[par].value = match[0].value
        return params

    measurements_loc = '../measurements.csv'

    loader = AngeliDataMapper(measurements_loc)
    analyzer = RadiusAnalyzer()
    grouper = DataGrouper()

    print('Absolute')
    rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'isoe', 'isom', 'noniso'])
    rad_info = analyzer.adjust_uncertainties(rad_info)
    opt_rad_info = RadiiInformation().join([analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    compare_radii_information('../outputs/absl_comparison.xlsx', [angeli_abs_info(), opt_rad_info],
                              prefixes=['angeli_abs', 'hunter_abs'])

    print('Noniso: R+')
    rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'isoe', 'isom'])
    rad_info = analyzer.adjust_uncertainties(rad_info)
    opt_rad_info = RadiiInformation().join(
        [analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    compare_radii_information('../outputs/noniso_comparison.xlsx', [angeli_remnoniso_info(), opt_rad_info],
                              prefixes=['angeli_noniso', 'hunter_noniso'])
    non_iso_info = opt_rad_info
    for data in non_iso_info.values(): # clear correlations, so values are are independent.
        data.correlations = defaultdict(dict)

    # print('REM: R++')
    # rad_info = loader.load_data(bad_types=['ois', 'kalpha'])
    # rad_info = analyzer.adjust_uncertainties(rad_info)
    # opt_rad_info = RadiiInformation().join(
    #     [analyzer.optimize_radii(rad_group, params=prev_params(rad_group, opt_rad_info)) for rad_group in grouper.group(rad_info)])
    #
    # compare_radii_information('../outputs/rem_comparison.xlsx', [angeli_rem_info(), opt_rad_info],
    #                           prefixes=['angeli_rem', 'hunter_rem'])

    print('REM: no noniso')
    rad_info = loader.load_data(bad_types=['ois', 'kalpha', 'noniso'])
    rad_info = analyzer.adjust_uncertainties(rad_info)
    rad_info = rad_info.join([non_iso_info])
    opt_rad_info = RadiiInformation().join(
        [analyzer.optimize_radii(rad_group, params=prev_params(rad_group, opt_rad_info)) for rad_group in
         grouper.group(rad_info)])

    compare_radii_information('../outputs/rem_without_nonisocomparison.xlsx', [angeli_rem_info(), opt_rad_info],
                              prefixes=['angeli_rem', 'hunter_rem'])