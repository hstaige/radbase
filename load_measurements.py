import pandas as pd

pd.set_option('display.max_columns', None)  # print everything
pd.set_option('display.width', 200)

iso_df = pd.DataFrame(columns=['Z', 'A'], dtype='int32')

abs_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', skiprows=5, usecols=list(range(1, 7)))
abs_df = abs_df.dropna(how='all')

cols = ['Z', 'A']
abs_df.loc[:, cols] = abs_df.loc[:, cols].ffill()
abs_df = abs_df[~abs_df['Z'].isin([0, 1])]
abs_df[['Type', 'Remarks']] = abs_df['Remarks'].str.extract(r'(el-scatt.|mu-Xray)\.*\s*(.*)', expand=True)
abs_df['Remarks'] = abs_df['Remarks'].fillna('')
abs_df = abs_df.astype({'Z': 'int32', 'A': 'int32', 'Ref.': 'string', 'Remarks': 'string'})
abs_df.rename(columns={'Rexp': 'R', 'dRexp': 'DR'}, inplace=True)
abs_df.reset_index(drop=True, inplace=True)
# abs_df.info(show_counts=True)

iso_df = pd.concat([iso_df, abs_df[['Z', 'A']]])

rel_noniso_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 5, 10, 13, 18, 31],
                              sheet_name='dRem_non_iso_11')
print(rel_noniso_df)
rel_noniso_df.drop(rel_noniso_df.tail(2).index, inplace=True)
rel_noniso_df.dropna(how='any', subset=['DRexp'], inplace=True)
cols = ['Z1', 'A1', 'Z2', 'A2']
rel_noniso_df.loc[:, cols] = rel_noniso_df.loc[:, cols].ffill()
rel_noniso_df.rename(columns={'DRexp': 'dR', 'd+D': 'DdR'}, inplace=True)
rel_noniso_df.reset_index(drop=True, inplace=True)
print('RELATIVE NONISOTOPIC')
rel_noniso_df.info(show_counts=True)

iso_df = pd.concat([iso_df, rel_noniso_df[['Z1', 'A1']].rename(columns={'Z1': 'Z', 'A1': 'A'}),
                    rel_noniso_df[['Z2', 'A2']].rename(columns={'Z2': 'Z', 'A2': 'A'})], axis=0)

rel_isoe_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 4, 7, 8, 26],
                            sheet_name='dRe_iso_123')
rel_isoe_df.dropna(how='any', subset=['DRexp'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_isoe_df.loc[:, cols] = rel_isoe_df.loc[:, cols].ffill()
rel_isoe_df.insert(3, 'Z2', rel_isoe_df['Z'])
rel_isoe_df.rename(columns={'Z': 'Z1', 'DRexp': 'dR', 'd+D': 'DdR'}, inplace=True)
print('\nRELATIVE ISOTOPIC ELECTRONIC')
rel_isoe_df.info(show_counts=True)
rel_isoe_df.reset_index(drop=True, inplace=True)

iso_df = pd.concat([iso_df, rel_isoe_df[['Z1', 'A1']].rename(columns={'Z1': 'Z', 'A1': 'A'}),
                    rel_isoe_df[['Z2', 'A2']].rename(columns={'Z2': 'Z', 'A2': 'A'})])


rel_isom_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[1, 2, 3, 12, 15, 20, 63],
                            sheet_name='dRm_iso_123')
rel_isom_df.dropna(how='any', subset=['DRx'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_isom_df.loc[:, cols] = rel_isom_df.loc[:, cols].ffill()
rel_isom_df.insert(3, 'Z2', rel_isom_df['Z'])
rel_isom_df.rename(columns={'Z': 'Z1', 'DRx': 'dR', 'd+D': 'DdR'}, inplace=True)
print('\nRELATIVE ISOTOPIC MUONIC')
rel_isom_df.info(show_counts=True)
rel_isom_df.reset_index(drop=True, inplace=True)

iso_df = pd.concat([iso_df, rel_isom_df[['Z1', 'A1']].rename(columns={'Z1': 'Z', 'A1': 'A'}),
                    rel_isom_df[['Z2', 'A2']].rename(columns={'Z2': 'Z', 'A2': 'A'})])

rel_k_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[1, 2, 3, 24, 29, 8, 56],
                         sheet_name='dRk_14')
rel_k_df.dropna(how='any', subset=['DR2'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_k_df.loc[:, cols] = rel_k_df.loc[:, cols].ffill()
rel_k_df.insert(3, 'Z2', rel_k_df['Z'])
rel_k_df.rename(columns={'Z': 'Z1', 'DR2': 'dR', 'dDR2': 'DdR'}, inplace=True)
print('\nRELATIVE K-Alpha')
rel_k_df.info(show_counts=True)
rel_k_df.reset_index(drop=True, inplace=True)

iso_df = pd.concat([iso_df, rel_k_df[['Z1', 'A1']].rename(columns={'Z1': 'Z', 'A1': 'A'}),
                    rel_k_df[['Z2', 'A2']].rename(columns={'Z2': 'Z', 'A2': 'A'})])

rel_ois_df = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 4, 7, 8, 9, 19],
                           sheet_name='dRo_19')
rel_ois_df.dropna(how='any', subset=['d<r2>'], inplace=True)
cols = ['Z', 'A1', 'A2']
rel_ois_df.loc[:, cols] = rel_ois_df.loc[:, cols].ffill()
rel_ois_df.insert(3, 'Z2', rel_ois_df['Z'])
rel_ois_df.rename(columns={'Z': 'Z1', 'd<r2>': 'dR', 'Dtot': 'DdR'}, inplace=True)
print('\nRELATIVE OIS')
rel_ois_df.info(show_counts=True)
rel_ois_df.reset_index(drop=True, inplace=True)

iso_df = pd.concat([iso_df, rel_ois_df[['Z1', 'A1']].rename(columns={'Z1': 'Z', 'A1': 'A'}),
                    rel_ois_df[['Z2', 'A2']].rename(columns={'Z2': 'Z', 'A2': 'A'})], axis=0)
iso_df.drop_duplicates(inplace=True)
iso_df.sort_values(by=['Z', 'A'], ascending=True, inplace=True)
iso_df.reset_index(drop=True, inplace=True)


def get_idx(iso, df, zcol, acol, idx_name='iso_idx'):
    """ Helper function that adds 'iso_idx' column to df"""
    idxs = [iso[(iso['Z'] == row[zcol]) & (iso['A'] == row[acol])].index[0] for (idx, row) in df.iterrows()]
    df.insert(0, idx_name, idxs)
    df.drop(columns=[zcol, acol], inplace=True)


get_idx(iso_df, abs_df, 'Z', 'A')
get_idx(iso_df, rel_noniso_df, 'Z2', 'A2', 'iso_idx2')
get_idx(iso_df, rel_noniso_df, 'Z1', 'A1', 'iso_idx1')
get_idx(iso_df, rel_isoe_df, 'Z2', 'A2', 'iso_idx2')
get_idx(iso_df, rel_isoe_df, 'Z1', 'A1', 'iso_idx1')
get_idx(iso_df, rel_isom_df, 'Z2', 'A2', 'iso_idx2')
get_idx(iso_df, rel_isom_df, 'Z1', 'A1', 'iso_idx1')
get_idx(iso_df, rel_k_df, 'Z2', 'A2', 'iso_idx2')
get_idx(iso_df, rel_k_df, 'Z1', 'A1', 'iso_idx1')
get_idx(iso_df, rel_ois_df, 'Z2', 'A2', 'iso_idx2')
get_idx(iso_df, rel_ois_df, 'Z1', 'A1', 'iso_idx1')


def add_rel_term(df, idx1, idx2, diff_type='linear'):
    if diff_type == 'linear':
        df['Term'] = 'R' + df[idx2].astype(str) + '-' + 'R' + df[idx1].astype(str)
    if diff_type == 'squared':
        df['Term'] = 'R' + df[idx2].astype(str) + '**2-' + 'R' + df[idx1].astype(str) + '**2'


abs_df['Term'] = 'R' + abs_df['iso_idx'].astype(str)
for df in [rel_noniso_df, rel_isoe_df, rel_isom_df, rel_k_df]:
    add_rel_term(df, 'iso_idx1', 'iso_idx2', diff_type='linear')
add_rel_term(rel_ois_df, 'iso_idx1', 'iso_idx2', diff_type='squared')


rel_noniso_df['Table'] = 'noniso'
rel_isoe_df['Table'] = 'isoe'
rel_isom_df['Table'] = 'isom'
rel_k_df['Table'] = 'kalpha'
rel_ois_df['Table'] = 'ois'


rel_df = pd.concat([rel_noniso_df, rel_isoe_df, rel_isom_df, rel_k_df, rel_ois_df], join='inner', ignore_index=True)
rel_df.to_csv('relative_radii.csv', index=False)
abs_df.to_csv('absolute_radii.csv', index=False)
iso_df.to_csv('isotope_indices.csv', index=False)