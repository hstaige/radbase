import pandas as pd
from radbase.solve import compress_measurements, save_human_readable


if __name__ == '__main__':

    pd.set_option('display.width', None)
    meas_df = pd.read_csv('../measurements.csv')

    meas_absonly, group_info_absonly = compress_measurements(measurements=meas_df[meas_df['Table'].isin(['absl'])])
    save_human_readable('../outputs/absonly', meas_absonly, group_info_absonly, 'absonly')

    meas_relonly, group_info_relonly = compress_measurements(
        measurements=meas_df[~meas_df['Table'].isin(['abs', 'noniso', 'ois'])])
    save_human_readable('../outputs/relonly', meas_relonly, group_info_relonly, 'relonly')

