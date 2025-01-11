import pandas as pd

from radbase.radius_analyzer import RadiiInformation


def get_idx(df, zcol, acol, idx_name='iso_idx'):
    """ Helper function that adds 'iso_idx' column to df"""
    idxs = 'R' + df[zcol].astype(int).astype(str).str.pad(width=3, side='left', fillchar='0') + df[acol].astype(int).astype(
        str).str.pad(width=3, side='left', fillchar='0')
    df.insert(0, idx_name, idxs)


def convert_to_radii_information(df):
    rad_info = RadiiInformation()
    for (i, row) in df.iterrows():
        rad_info.add(term=row['Term'], value=row['Value'], unc=row['Unc'])
    return rad_info


def angeli_abs_info(file_loc='../NuclearRadius_all_sheets_renamed.xlsx'):
    abs_df = pd.read_excel(file_loc, skiprows=5)
    abs_df[['Z', 'A']] = abs_df[['Z', 'A']].ffill()
    abs_df = abs_df.dropna(subset='Rav')
    get_idx(abs_df, 'Z', 'A', idx_name='Term')
    abs_df.rename(columns={'Rav': 'Value', 'dRav': 'Unc'}, inplace=True)
    return convert_to_radii_information(abs_df)


def angeli_remnoniso_info(file_loc='../NuclearRadius_all_sheets_renamed.xlsx'):
    abs_df = pd.read_excel(file_loc, sheet_name='Rem_non_iso_11')
    abs_df = abs_df.dropna(subset='R+')
    abs_df[['Z', 'A']] = abs_df[['Z', 'A']].ffill()
    abs_df = abs_df[abs_df['Include'] != 'No']
    get_idx(abs_df, 'Z', 'A', idx_name='Term')
    abs_df.rename(columns={'R+': 'Value', 'DR+': 'Unc'}, inplace=True)
    abs_rad_info = convert_to_radii_information(abs_df)

    rel_df = pd.read_excel(file_loc, sheet_name='Rem_non_iso_11')
    rel_df = rel_df[rel_df['Include'] != 'No']
    rel_df = rel_df.dropna(subset='dR+')
    get_idx(rel_df, 'Z1', 'A1', idx_name='iso1')
    get_idx(rel_df, 'Z2', 'A2', idx_name='iso2')
    rel_df['Term'] = rel_df['iso2'] + '-' + rel_df['iso1']
    rel_df.rename(columns={'dR+': 'Value', 'DdR+': 'Unc'}, inplace=True)
    rel_rad_info = convert_to_radii_information(rel_df)

    return RadiiInformation().join([rel_rad_info, abs_rad_info])


def angeli_rem_info(file_loc='../NuclearRadius_all_sheets_renamed.xlsx'):
    abs_df = pd.read_excel(file_loc, sheet_name='Rem_151')
    abs_df = abs_df.dropna(subset='R++')
    abs_df.loc[abs_df['Z'].astype(str).str.contains('!'), 'Z'] = pd.NA
    abs_df[['Z', 'A']] = abs_df[['Z', 'A']].ffill()
    abs_df = abs_df[abs_df['Include'] != 'No']
    get_idx(abs_df, 'Z', 'A', idx_name='Term')
    abs_df.rename(columns={'R++': 'Value', 'DR++': 'Unc'}, inplace=True)
    abs_rad_info = convert_to_radii_information(abs_df)

    rel_df = pd.read_excel(file_loc, sheet_name='Rem_151')
    rel_df.loc[rel_df['Z'].astype(str).str.contains('!'), 'Z'] = pd.NA
    rel_df[['Z', 'A']] = rel_df[['Z', 'A']].ffill()
    rel_df = rel_df.dropna(subset='dRem')
    rel_df = rel_df[rel_df['Include'] != 'No']
    get_idx(rel_df, 'Z', 'A1', idx_name='iso1')
    get_idx(rel_df, 'Z', 'A2', idx_name='iso2')
    rel_df['Term'] = rel_df['iso2'] + '-' + rel_df['iso1']
    rel_df.rename(columns={'dR++': 'Value', 'DdR++': 'Unc'}, inplace=True)
    rel_rad_info = convert_to_radii_information(rel_df)

    return RadiiInformation().join([rel_rad_info, abs_rad_info])

def angeli_remko_info(file_loc='../NuclearRadius_all_sheets_renamed.xlsx'):
    abs_df = pd.read_excel(file_loc, sheet_name='Remko_152')
    abs_df = abs_df.dropna(subset='R++')
    abs_df.loc[abs_df['Z'].astype(str).str.contains('!'), 'Z'] = pd.NA
    abs_df[['Z', 'A']] = abs_df[['Z', 'A']].ffill()
    abs_df = abs_df[abs_df['Include'] != 'No']
    get_idx(abs_df, 'Z', 'A', idx_name='Term')
    abs_df.rename(columns={'R++': 'Value', 'DR++': 'Unc'}, inplace=True)
    abs_rad_info = convert_to_radii_information(abs_df)

    rel_df = pd.read_excel(file_loc, sheet_name='Remko_152')
    rel_df.loc[rel_df['Z'].astype(str).str.contains('!'), 'Z'] = pd.NA
    rel_df[['Z', 'A']] = rel_df[['Z', 'A']].ffill()
    rel_df = rel_df.dropna(subset='dR++')
    rel_df = rel_df[rel_df['Include'] != 'No']
    get_idx(rel_df, 'Z', 'A1', idx_name='iso1')
    get_idx(rel_df, 'Z', 'A2', idx_name='iso2')
    rel_df['Term'] = rel_df['iso2'] + '-' + rel_df['iso1']
    rel_df.rename(columns={'dR++': 'Value', 'DdR++': 'Unc'}, inplace=True)
    rel_rad_info = convert_to_radii_information(rel_df)

    return RadiiInformation().join([rel_rad_info, abs_rad_info])


def angeli_rfinal_info(file_loc='../NuclearRadius_all_sheets_renamed.xlsx'):
    abs_df = pd.read_excel(file_loc, sheet_name='R_final_152')
    abs_df = abs_df.dropna(subset='R')
    abs_df[['Z', 'A']] = abs_df[['Z', 'A']].ffill()
    get_idx(abs_df, 'Z', 'A', idx_name='Term')
    abs_df.rename(columns={'R': 'Value', 'DtotR': 'Unc'}, inplace=True)
    abs_rad_info = convert_to_radii_information(abs_df)
    return abs_rad_info


if __name__ == '__main__':
    file_loc = '../../inputs/NuclearRadius_all_sheets_renamed.xlsx'
    # print(angeli_abs_info(file_loc))
    # print(angeli_remnoniso_info(file_loc))
    print(angeli_final_info(file_loc)[155])
    pass
