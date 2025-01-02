from radbase.radius_analyzer import RadiiInformation, RadiusAnalyzer, DataGrouper
from radbase.database import AngeliDataMapper


if __name__ == '__main__':

    measurements_loc = '../measurements.csv'

    loader = AngeliDataMapper(measurements_loc)
    analyzer = RadiusAnalyzer()
    grouper = DataGrouper()

    rad_info = loader.load_data(bad_types=['ois', 'kalpha'])
    opt_rad_info = RadiiInformation().join([analyzer.optimize_radii(rad_group) for rad_group in grouper.group(rad_info)])
    print(opt_rad_info)
