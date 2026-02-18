import pandas as pd

dfa = pd.read_csv('./cross_sections_a.csv')
dfb = pd.read_csv('./cross_sections_b.csv')
dfc = pd.read_csv('./cross_sections_c.csv')

default_angle_acceptance = 1.8
angle_acceptances = {'ᵃ': 0.9, 'ᵇ': 0.3, 'ᶜ': 0.36, 'ᵈ': 0.2}

if 'Energy [keV]' not in dfb.columns:
    dfb.insert(loc=0, column='Energy [keV]', value=250)
    dfb.to_csv('./cross_sections_b.csv', index=False)

for_all_cols = ['Energy [keV]', 'theta']

by_quant_dfs = []
for df in [dfa, dfb, dfc]:
    for column in df.columns:
        if column not in for_all_cols:
            col_df = df[~df[column].isin(['', '...'])]
            col_df = col_df[for_all_cols + [column]]
            col_df.rename(columns={column: 'Value'}, inplace=True)
            col_df['Quantity'] = column
            by_quant_dfs.append(col_df)

by_quant_dfs = pd.concat(by_quant_dfs)
by_quant_dfs.to_csv('./cross_sections.csv', index=False)

for i, row in by_quant_dfs.iterrows():
    print(row)
