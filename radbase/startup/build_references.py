import json
import os

import pandas as pd
from crossref_commons.iteration import iterate_publications_as_json
from crossref_commons.retrieval import get_publication_as_json

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


def query_crossref(citation, num_return=3):
    # iter_json things
    filter = {'type': 'journal-article'}
    queries = {'query.bibliographic': citation}
    return list(iterate_publications_as_json(max_results=num_return, filter=filter, queries=queries))


def query_driver(df, num_per_backup=1, **kwargs):

    def display_json(i, opt):
        try:
            print(f'{i:4}: Title: {opt['title'][0]}\n'
                  f'      Published {opt['published']['date-parts'][0]} in {opt['short-container-title'][0]}\n'
                  f'      By: {', '.join([auth.get('given', '') + ' ' + auth.get('family', '') for auth in opt.get('author', [{'given': 'missing'}])])}\n')
        except KeyError as e:
            print(f'{i}: missing {e}')

    def save_ref_json(ref, refid):
        path = f'.//references//{refid}.json'
        if os.path.isfile(path):
            answer = input('This file already exists! Overwrite? (y/n) ')
            if answer != 'y':
                return

        with open(path, 'w') as f:
            json.dump(ref, f)

    skipto = int(input(f'Which index out of {len(df)} would you like to begin at? '))

    selected = {}
    for row_idx, row in df[skipto:].iterrows():
        print(f'Querying on {row_idx}/{len(df)}, {row['id']} - {row['Citation']} ... \n')
        options = query_crossref(row['Citation'], **kwargs)
        for i, option in enumerate(options):
            display_json(i, option)
        selection = input('Select an index from the above, -1 to skip, or doi to enter doi: ')

        if selection == 'doi':
            doi = input('Please enter doi! ')
            try:
                selected[row['id']] = get_publication_as_json(doi)
            except:
                print('DOI retrieval failed ...')
        else:
            try:
                selection = int(selection)
            except ValueError:
                print('Invalid index?')
                selection = -1
            if selection == -1:
                print('Skipping ... \n')
                continue

            selected[row['id']] = options[selection]

        if len(selected) >= num_per_backup:
            for ref_id in selected:
                save_ref_json(selected[ref_id], ref_id)
            selected = {}




# ref_df = pd.read_csv(r'../../references.csv')
# query_driver(ref_df)


# res = process_garbage('../../References.xlsx')
# print(res.info())
# res.to_csv('../../references.csv', index=False)
