"""
Converts Angeli's Excel file into the required input format described in docs/inputs.rst

Last updated 12/14/24

"""

import re

import pandas as pd

pd.set_option('display.max_columns', None)  # print everything
pd.set_option('display.width', 200)

abs_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', skiprows=5, usecols=list(range(1, 7)))
abs_df = abs_df.dropna(how='all')

cols = ['Z', 'A']
abs_df.loc[:, cols] = abs_df.loc[:, cols].ffill()
abs_df = abs_df[~abs_df['Z'].isin([0, 1])]
abs_df = abs_df.astype({'Z': 'int32', 'A': 'int32', 'Ref.': 'string', 'Remarks': 'string'})
abs_df['RefRemarks'] = abs_df['Ref.'] + ' ' + abs_df['Remarks']
abs_df['Refs'] = abs_df['RefRemarks'].apply(lambda r: re.findall(r'\w\w\d\d[a-e]?', r))
print(abs_df['Refs'][abs_df['Refs'].notna()])
abs_df.drop(columns=['Ref.', 'Remarks'], inplace=True)
abs_df.rename(columns={'Rexp': 'Value', 'dRexp': 'Unc'}, inplace=True)
abs_df.reset_index(drop=True, inplace=True)
# abs_df.info(show_counts=True)

rel_noniso_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 5, 10, 13, 18, 31],
                              sheet_name='dRem_non_iso_11')
rel_noniso_df.drop(rel_noniso_df.tail(2).index, inplace=True)
rel_noniso_df.dropna(how='any', subset=['DRexp'], inplace=True)
cols = ['Z1', 'A1', 'Z2', 'A2']
rel_noniso_df.loc[:, cols] = rel_noniso_df.loc[:, cols].ffill()
rel_noniso_df.rename(columns={'DRexp': 'Value', 'd+D': 'Unc', 'Ref. + Remarks': 'RefRemarks'}, inplace=True)
rel_noniso_df['Refs'] = rel_noniso_df['RefRemarks'].apply(lambda r: re.findall(r'(?<=\()[a-zA-Z]{2}\d\d[a-e]?', r) +
                                                                    re.findall(r'(?<!\()[a-zA-Z]{2}\d\d[a-e]?', r))
print(rel_noniso_df['Refs'])
rel_noniso_df.reset_index(drop=True, inplace=True)
print('RELATIVE NONISOTOPIC')
rel_noniso_df.info(show_counts=True)

rel_isoe_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 4, 7, 8, 26],
                            sheet_name='dRe_iso_123')
rel_isoe_df.dropna(how='any', subset=['DRexp'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_isoe_df.loc[:, cols] = rel_isoe_df.loc[:, cols].ffill()
rel_isoe_df.insert(3, 'Z2', rel_isoe_df['Z'])
rel_isoe_df.rename(columns={'Z': 'Z1', 'DRexp': 'Value', 'd+D': 'Unc', 'Ref.': 'RefRemarks'}, inplace=True)
rel_isoe_df['Refs'] = rel_isoe_df['RefRemarks'].apply(lambda r: re.findall(r'(?<=\()[a-zA-Z]{2}\d\d[a-e]?', r) +
                                                                re.findall(r'(?<!\()[a-zA-Z]{2}\d\d[a-e]?', r))
print('\nRELATIVE ISOTOPIC ELECTRONIC')
rel_isoe_df.info(show_counts=True)
rel_isoe_df.reset_index(drop=True, inplace=True)

rel_isom_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', usecols=[1, 2, 3, 12, 15, 20, 63],
                            sheet_name='dRm_iso_123')
rel_isom_df.dropna(how='any', subset=['DRx'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_isom_df.loc[:, cols] = rel_isom_df.loc[:, cols].ffill()
rel_isom_df.insert(3, 'Z2', rel_isom_df['Z'])
rel_isom_df.rename(columns={'Z': 'Z1', 'DRx': 'Value', 'd+D': 'Unc', 'Ref. + Rem.': 'RefRemarks'}, inplace=True)
rel_isom_df['Refs'] = rel_isom_df['RefRemarks'].apply(lambda r: re.findall(r'(?<=\()[a-zA-Z]{2}\d\d[a-e]?', r) +
                                                                re.findall(r'(?<!\()[a-zA-Z]{2}\d\d[a-e]?', r))
print('\nRELATIVE ISOTOPIC MUONIC')
rel_isom_df.info(show_counts=True)
rel_isom_df.reset_index(drop=True, inplace=True)

rel_k_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', usecols=[1, 2, 3, 24, 29, 8, 56],
                         sheet_name='dRk_14')
rel_k_df.dropna(how='any', subset=['DR2'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_k_df.loc[:, cols] = rel_k_df.loc[:, cols].ffill()
rel_k_df.insert(3, 'Z2', rel_k_df['Z'])
rel_k_df.rename(columns={'Z': 'Z1', 'DR2': 'Value', 'dDR2': 'Unc', 'Ref.': 'RefRemarks'}, inplace=True)
rel_k_df['Refs'] = rel_k_df['RefRemarks'].apply(lambda r: re.findall(r'(?<=\()[a-zA-Z]{2}\d\d[a-e]?', r) +
                                                                re.findall(r'(?<!\()[a-zA-Z]{2}\d\d[a-e]?', r))
print('\nRELATIVE K-Alpha')
rel_k_df.info(show_counts=True)
rel_k_df.reset_index(drop=True, inplace=True)

rel_ois_df = pd.read_excel('../../inputs/NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 4, 7, 8, 9, 19],
                           sheet_name='dRo_19')
rel_ois_df.dropna(how='any', subset=['d<r2>'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_ois_df.loc[:, cols] = rel_ois_df.loc[:, cols].ffill()
rel_ois_df.insert(3, 'Z2', rel_ois_df['Z'])
rel_ois_df['RefRemarks'] = rel_ois_df['Ref.'].astype(str) + rel_ois_df['Remark'].astype(str)
rel_ois_df.rename(columns={'Z': 'Z1', 'd<r2>': 'Value', 'Dtot': 'Unc'}, inplace=True)
rel_ois_df['Refs'] = None
print('\nRELATIVE OIS')
rel_ois_df.info(show_counts=True)
rel_ois_df.reset_index(drop=True, inplace=True)


def get_idx(df, zcol, acol, idx_name='iso_idx'):
    """ Helper function that adds 'iso_idx' column to df"""
    idxs = df[zcol].astype(int).astype(str).str.pad(width=3, side='left', fillchar='0') + df[acol].astype(int).astype(
        str).str.pad(width=3, side='left', fillchar='0')
    df.insert(0, idx_name, idxs)
    df.drop(columns=[zcol, acol], inplace=True)


get_idx(abs_df, 'Z', 'A')
get_idx(rel_noniso_df, 'Z2', 'A2', 'iso_idx2')
get_idx(rel_noniso_df, 'Z1', 'A1', 'iso_idx1')
get_idx(rel_isoe_df, 'Z2', 'A2', 'iso_idx2')
get_idx(rel_isoe_df, 'Z1', 'A1', 'iso_idx1')
get_idx(rel_isom_df, 'Z2', 'A2', 'iso_idx2')
get_idx(rel_isom_df, 'Z1', 'A1', 'iso_idx1')
get_idx(rel_k_df, 'Z2', 'A2', 'iso_idx2')
get_idx(rel_k_df, 'Z1', 'A1', 'iso_idx1')
get_idx(rel_ois_df, 'Z2', 'A2', 'iso_idx2')
get_idx(rel_ois_df, 'Z1', 'A1', 'iso_idx1')


def add_rel_term(df, idx1, idx2, diff_type='linear'):
    if diff_type == 'linear':
        df['Term'] = 'R' + df[idx2].astype(str) + '-' + 'R' + df[idx1].astype(str)
    if diff_type == 'squared':
        df['Term'] = 'R' + df[idx2].astype(str) + '**2-' + 'R' + df[idx1].astype(str) + '**2'


abs_df['Term'] = 'R' + abs_df['iso_idx'].astype(str)
for df in [rel_noniso_df, rel_isoe_df, rel_isom_df]:
    add_rel_term(df, 'iso_idx1', 'iso_idx2', diff_type='linear')
add_rel_term(rel_k_df, 'iso_idx1', 'iso_idx2', diff_type='squared')
add_rel_term(rel_ois_df, 'iso_idx1', 'iso_idx2', diff_type='squared')

rel_noniso_df['Table'] = 'noniso'
rel_isoe_df['Table'] = 'isoe'
rel_isom_df['Table'] = 'isom'
rel_k_df['Table'] = 'kalpha'
rel_ois_df['Table'] = 'ois'

abs_df.rename(columns={'iso_idx': 'Iso_Idxs'}, inplace=True)
abs_df = abs_df.astype({'Iso_Idxs': str})
abs_df['Include'] = ''
abs_df['Table'] = 'absl'

rel_df = pd.concat([rel_noniso_df, rel_isoe_df, rel_isom_df, rel_k_df, rel_ois_df], join='inner', ignore_index=True)
rel_df['Iso_Idxs'] = rel_df['iso_idx1'] + ',' + rel_df['iso_idx2']
print(rel_df)
rel_df.drop(columns=['iso_idx1', 'iso_idx2'], inplace=True)

measurement_df = pd.concat([abs_df, rel_df], ignore_index=True)
measurement_df.reset_index()

measurement_df = measurement_df.astype({'Iso_Idxs': str, 'Value': float, 'Unc': float,
                                        'RefRemarks': str, 'Term': str, 'Include': str,
                                        'Table': str})

measurement_df.info(show_counts=True)

measurement_df.to_csv('../../inputs/measurements.csv', index=False)
