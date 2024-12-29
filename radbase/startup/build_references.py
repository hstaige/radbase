import pandas as pd

pd.set_option('display.max_columns', None)  # print everything
pd.set_option('display.width', None)

def process_garbage(file_loc):
    df = pd.read_excel(file_loc)
    df[['id', 'Description']] = df['Description'].str.split(' ', n=1, expand=True)
    df[['Citation', 'Elements']] = df['Description'].str.split('..', n=1, regex=False, expand=True)
    df.drop(columns=['Description'], inplace=True)
    df['Authors'] = df['Citation'].str.split(':', n=1, regex=False, expand=True)[0]
    df['Authors'] = df['Authors'].str.replace(' and ', '').str.replace(' ', '')
    df['Year'] = df['Citation'].str.extract(r'(\d{4})')

    return df


res = process_garbage('../../References.xlsx')
print(res.info())
res.to_csv('references.csv', index=False)