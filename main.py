import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import lmfit
import asteval
from itertools import chain

pd.set_option('display.width', 300)

rel_df = pd.read_csv('relative_radii.csv')
rel_df = rel_df[rel_df['Include'] != 'No']
abs_df = pd.read_csv('absolute_radii.csv')
iso_df = pd.read_csv('isotope_indices.csv')


def compress_measurements(rel=None, absl=None):
    # idea: start with absolute measurements, build out connected groups using rel to connect
    #       any left over relative measurements should be connected, treated like floating?

    # therefore we want bfs floodfill that removes edges as they are crossed
    # each edge should store set of measurements corresponding to that edge?

    # graph[idx] -> {nextidx: {rel_idxs}, ...}

    def update_inputs(grp, sol, r, a):

        rel_idxs = r.loc[grp['rel_idxs'], 'connect_rel_idxs'].to_list() + a.loc[grp['abs_idxs'], 'connect_rel_idxs'].to_list()
        rel_idxs = filter(None, rel_idxs)  # remove any blank strings
        rel_idxs = ','.join(set(','.join(rel_idxs).split(',')))  # combine list of comma separated idxs into one str
        abs_idxs = r.loc[grp['rel_idxs'], 'connect_abs_idxs'].to_list() + a.loc[grp['abs_idxs'], 'connect_abs_idxs'].to_list()
        abs_idxs = filter(None, abs_idxs)  # remove any blank strings
        abs_idxs = ','.join(set(','.join(abs_idxs).split(',')))  # combine list of comma separated idxs into one str

        multiplier = max(1, (0.5 * (1 + sol.redchi)) ** (1/2) if sol.redchi else 0)

        a.loc[grp['abs_idxs'], 'R'] = [sol.params[f'R{a.loc[aidx, "iso_idx"]}'].value for aidx in grp['abs_idxs']]
        a.loc[grp['abs_idxs'], 'DR'] = [sol.params[f'R{a.loc[aidx, "iso_idx"]}'].stderr * multiplier for aidx in grp['abs_idxs']]
        a.loc[grp['abs_idxs'], 'connect_abs_idxs'] = abs_idxs
        a.loc[grp['abs_idxs'], 'connect_rel_idxs'] = rel_idxs

        r.loc[grp['rel_idxs'], 'dR'] = [sol.params[f'rel{ridx}'].value for ridx in grp['rel_idxs']]
        r.loc[grp['rel_idxs'], 'DdR'] = [sol.params[f'rel{ridx}'].stderr * multiplier for ridx in grp['rel_idxs']]
        r.loc[grp['rel_idxs'], 'connect_abs_idxs'] = abs_idxs
        r.loc[grp['rel_idxs'], 'connect_rel_idxs'] = rel_idxs
        return r, a

    def update_group_info(grp, sol, ginfo):
        # expected columns of ginfo: abs_idxs, rel_idxs, iso_idxs dof, redchisq, lim
        grp = {key: ','.join(map(str, list(el)) )for key, el in grp.items()}
        grp.update({'redchisq': sol.redchi, 'dof': sol.nfree, 'lim': 1 + 2 * (2 / sol.nfree) ** (1 / 2) if sol.nfree else 100})
        return pd.concat([ginfo, pd.DataFrame(grp, index=[len(ginfo)])])

    if rel is None and absl is None:
        raise ValueError('At least one of rel and absl must not be None')

    rel = rel if rel is not None else pd.DataFrame(columns=['iso_idx1', 'iso_idx2', 'dR', 'DdR', 'Term'])
    absl = absl if absl is not None else pd.DataFrame(columns=['iso_idx', 'R', 'DR', 'Term'])

    groups = find_groups(rel=rel, absl=absl)  # break up relative and absolute matrices into solid groups
    solutions = {i: solve_group(group, rel=rel, absl=absl) for i, group in groups.items()}

    absl, rel = absl.copy(), rel.copy()

    if 'connect_abs_idxs' not in absl.columns:
        absl['connect_abs_idxs'] = [str(i) for i in absl.index.to_list()]
    if 'connect_rel_idxs' not in absl.columns:
        absl['connect_rel_idxs'] = ['' for _ in absl.index.to_list()]
    if 'connect_abs_idxs' not in rel.columns:
        rel['connect_abs_idxs'] = ['' for _ in rel.index.to_list()]
    if 'connect_rel_idxs' not in rel.columns:
        rel['connect_rel_idxs'] = [str(i) for i in rel.index.to_list()]

    grp_info = pd.DataFrame(columns=['abs_idxs', 'rel_idxs', 'iso_idxs', 'dof', 'redchisq', 'lim'])
    for group, solution in zip(groups.values(), solutions.values()):
        rel, absl = update_inputs(group, solution, rel, absl)
        grp_info = update_group_info(group, solution, grp_info)

    rel.drop_duplicates(subset=['iso_idx1', 'iso_idx2'], inplace=True)
    absl.drop_duplicates(subset=['iso_idx'], inplace=True)

    return rel, absl, grp_info


def find_groups(rel=None, absl=None):
    if rel is None and absl is None:
        raise ValueError('At least one of rel and absl must not be None')

    graph = defaultdict(dict)
    adict = defaultdict(set)
    if absl is not None:
        for iso_idx, aidx in zip(absl['iso_idx'], absl.index):
            adict[iso_idx].add(aidx)

    if rel is not None:
        # graph[idx1] -> {idx2: {rel_idxs}, ...}
        for r1, r2, ridx in zip(rel['iso_idx1'], rel['iso_idx2'], rel.index):
            if r2 not in graph[r1]:  # we haven't hit this edge yet, add it to the graph
                graph[r1][r2] = set()
                graph[r2][r1] = set()

            graph[r1][r2].add(ridx)
            graph[r2][r1].add(ridx)

    rel_idxs = [] if rel is None else rel['iso_idx1'].tolist() + rel['iso_idx2'].to_list()
    abs_idxs = [] if absl is None else absl['iso_idx'].to_list()

    seen = set()
    groups = []
    isotopes = set(rel_idxs + abs_idxs)
    for i in isotopes:
        if i in seen:
            continue

        group = {'abs_idxs': set(), 'rel_idxs': set(), 'iso_idxs': set()}
        stack = [i]
        seen.add(i)
        while stack:
            curr = stack.pop()
            group['abs_idxs'].update(adict[curr])  # add all absolute measurements for this isotope
            group['iso_idxs'].add(curr)

            for next_i in list(graph[curr].keys()):  # cant just do graph[curr] because we want to delete elements
                group['rel_idxs'].update(graph[curr][next_i])  # add rel_idxs for this edge to the group
                del graph[curr][next_i]  # delete the edge of the graph
                if next_i in seen:
                    continue
                seen.add(next_i)  # mark as already visited
                stack.append(next_i)

        groups.append(group)

    return {i: group for i, group in enumerate(groups)}


def solve_group(group, rel, absl):
    params = lmfit.Parameters()

    for iso in group['iso_idxs']:
        params.add(f'R{iso}', 0)

    for iso, R in absl.loc[group['abs_idxs'], ['iso_idx', 'R']].itertuples(index=False):
        params[f'R{iso}'].value = R

    for ridx, term in zip(group['rel_idxs'], rel.loc[group['rel_idxs'], 'Term']):
        params.add(f'rel{ridx}', expr=term)

    targets = list(absl.loc[group['abs_idxs'], ['Term', 'R', 'DR']].itertuples(index=False))
    targets += list(rel.loc[group['rel_idxs'], ['Term', 'dR', 'DdR']].itertuples(index=False))

    nfree = sum(params[p].vary for p in params)
    terms = list(zip(*targets))[0]
    if not group['abs_idxs'] and (all('**2' not in t for t in terms) or (nfree > len(targets))):  # no absolute measurements and only linear measurements
        params[list(params.keys())[0]].vary = False  # fix one of the isotopes so variance is not infinite

    def target(pars, tars):
        return [(pars.eval(term) - val) / unc for term, val, unc in tars]

    mini_res = lmfit.minimize(target, params, args=(targets,), scale_covar=False)

    if 20 in group['iso_idxs']:
        print(mini_res.params)
        print(targets)

    return mini_res


def make_human_readable(rel, absl, iso):
    rel, absl = rel.copy(), absl.copy()

    print(rel.columns)
    if rel.columns[0] != 'iso_idx1':
        raise Warning('relative columns may be out of order')
    rel.insert(0, 'Z1', iso.loc[rel['iso_idx1'], 'Z'].to_list())
    rel.insert(1, 'A1', iso.loc[rel['iso_idx1'], 'A'].to_list())
    rel.insert(2, 'Z2', iso.loc[rel['iso_idx2'], 'Z'].to_list())
    rel.insert(3, 'A2', iso.loc[rel['iso_idx2'], 'A'].to_list())
    rel.drop(columns=['iso_idx1', 'iso_idx2'], inplace=True)
    rel.sort_values(by=['Z1', 'A1', 'Z2', 'A2'], inplace=True)
    rel = rel.astype({'Z1': int, 'A1': int, 'Z2': int, 'A2': int})

    if absl.columns[0] != 'iso_idx':
        raise Warning('absolute columns may be out of order')
    absl.insert(0, 'Z', iso.loc[absl['iso_idx'], 'Z'].to_list())
    absl.insert(2, 'A', iso.loc[absl['iso_idx'], 'A'].to_list())
    absl.drop(columns=['iso_idx', 'Term'], inplace=True)
    absl.sort_values(by=['Z', 'A'], inplace=True)
    absl = absl.astype({'Z': int, 'A': int})

    return rel, absl


# print(find_groups(rel=rel_df))
# print(find_groups(absl=abs_df))
# print(find_groups(rel=rel_df, absl=abs_df))

rel_absonly, abs_absonly, group_info_absonly = compress_measurements(rel=None, absl=abs_df)
rel_absonly_hr, abs_absonly_hr = make_human_readable(rel_absonly, abs_absonly, iso_df)
abs_absonly_hr.to_csv('./absonly/absolute_radii_absonly.csv', index=False)
rel_absonly_hr.to_csv('./absonly/relative_radii_absonly.csv', index=False)
group_info_absonly.to_csv('./absonly/group_info_absonly.csv', index=False)

rel_relonly, abs_relonly, group_info_relonly = compress_measurements(rel=rel_df[~rel_df['Table'].isin(['noniso', 'ois'])], absl=None)
rel_relonly, abs_relonly = make_human_readable(rel_relonly, abs_relonly, iso_df)
abs_relonly.to_csv('./relonly/absolute_radii_relonly.csv', index=False)
rel_relonly.to_csv('./relonly/relative_radii_relonly.csv', index=False)
group_info_relonly.to_csv('./relonly/group_info_relonly.csv', index=False)

iso_copy = iso_df.copy().reset_index() # make index of iso_df its own column named #index
rel_slice = pd.read_excel('./NuclearRadius_all_sheets_renamed.xlsx', usecols=[0, 2, 3, 5, 32],
                              sheet_name='dRem_non_iso_11')
rel_slice = rel_slice[rel_slice['AdjUnc'] == 'Yes']
rel_slice = pd.merge(iso_copy, rel_slice, left_on=['Z', 'A'], right_on=['Z1', 'A1'])
rel_slice.rename(columns={'index': 'iso_idx1'}, inplace=True)
rel_slice = pd.merge(iso_copy, rel_slice, left_on=['Z', 'A'], right_on=['Z2', 'A2'])
rel_slice.rename(columns={'index': 'iso_idx2'}, inplace=True)


rel_noniso_relonly, abs_noniso_relonly, group_info_noniso_relonly = compress_measurements(rel=rel_df[rel_df['Table'] == 'noniso'], absl=None)
rel_noniso_relonly_hr, abs_noniso_relonly_hr = make_human_readable(rel_noniso_relonly, abs_noniso_relonly, iso_df)
abs_noniso_relonly_hr.to_csv('./noniso_relonly/absolute_radii_noniso_relonly.csv', index=False)
rel_noniso_relonly_hr.to_csv('./noniso_relonly/relative_radii_noniso_relonly.csv', index=False)
group_info_noniso_relonly.to_csv('./noniso_relonly/group_info_noniso_relonly.csv', index=False)

rel2 = rel_noniso_relonly.copy()
rel1 = rel_df[rel_df['Table'] == 'noniso']
for i, row in rel_slice.iterrows():
    i1, i2 = row['iso_idx1'], row['iso_idx2']
    orig_row = rel1[(rel1['iso_idx1'] == i1) & (rel1['iso_idx2'] == i2)].iloc[0]
    new_idx = rel_noniso_relonly[(rel_noniso_relonly['iso_idx1'] == i1) & (rel_noniso_relonly['iso_idx2'] == i2)].index[0]
    rel_noniso_relonly.loc[new_idx, 'DdR'] = orig_row['DdR']

print(rel_noniso_relonly.head(20))

rel_noniso, abs_noniso, group_info_noniso = compress_measurements(rel=rel_noniso_relonly, absl=abs_absonly)
rel_noniso_hr, abs_noniso_hr = make_human_readable(rel_noniso, abs_noniso, iso_df)
abs_noniso_hr.to_csv('./noniso/absolute_radii_noniso.csv', index=False)
rel_noniso_hr.to_csv('./noniso/relative_radii_noniso.csv', index=False)
group_info_noniso.to_csv('./noniso/group_info_noniso.csv', index=False)

print(rel_df[~rel_df['Table'].isin(['noniso', 'ois', 'kalpha'])]['Table'].unique())

# rel_rem, abs_rem, group_info_rem = compress_measurements(rel=rel_df[~rel_df['Table'].isin(['noniso', 'ois', 'kalpha'])], absl=abs_noniso)
# rel_rem, abs_rem = make_human_readable(rel_rem, abs_rem, iso_df)
# abs_rem.to_csv('./rem/absolute_radii_rem.csv', index=False)
# rel_rem.to_csv('./rem/relative_radii_rem.csv', index=False)
# group_info_rem.to_csv('./rem/group_info_rem.csv', index=False)
#
#
# rel_all, abs_all, group_info_all = compress_measurements(rel=rel_df[~(rel_df['Table'] == 'noniso')], absl=abs_noniso)
#
# rel_all, abs_all = make_human_readable(rel_all, abs_all, iso_df)
# abs_all.to_csv('./all/absolute_radii_all.csv', index=False)
# rel_all.to_csv('./all/relative_radii_all.csv', index=False)
# group_info_all.to_csv('./all/group_info_all.csv', index=False)
