import pandas as pd
from collections import defaultdict
from itertools import product, chain
import matplotlib.pyplot as plt
import lmfit
import asteval
from itertools import chain


def compress_measurements(measurements):
    def update_inputs(grp, sol, m):
        midxs = grp['m_idxs']

        connected = set([a for b in m.loc[midxs, 'Connected_Measurements'].str.split(',').tolist() for a in b])

        m.loc[midxs, 'Connected_Measurements'] = ','.join(list(connected))

        multiplier = max(1, (0.5 * (1 + sol.redchi)) ** (1 / 2) if sol.redchi else 0)
        m.loc[midxs, 'Value'] = [sol.params[f'm{midx}'].value for midx in midxs]
        m.loc[midxs, 'Unc'] = [sol.params[f'm{midx}'].stderr * multiplier for midx in midxs]

        return m

    def update_group_info(grp, sol, ginfo):
        # expected columns of ginfo: abs_idxs, rel_idxs, iso_idxs dof, redchisq, lim
        grp = {key: ','.join(map(str, list(el))) for key, el in grp.items()}
        grp.update(
            {'redchisq': sol.redchi, 'dof': sol.nfree, 'lim': 1 + 2 * (2 / sol.nfree) ** (1 / 2) if sol.nfree else 100})
        return pd.concat([ginfo, pd.DataFrame(grp, index=[len(ginfo)])])

    measurements = measurements.copy()

    groups = find_groups(measurements)
    solutions = [solve_group(group, measurements) for group in groups]

    if 'Connected_Measurements' not in measurements.columns:
        measurements['Connected_Measurements'] = measurements.index
        measurements = measurements.astype({'Connected_Measurements': str})

    grp_info = pd.DataFrame()
    for idx, (group, solution) in enumerate(zip(groups, solutions)):
        measurements = update_inputs(group, solution, measurements)
        grp_info = update_group_info(group, solution, grp_info)

    measurements.drop_duplicates('Iso_Idxs')

    return measurements, grp_info


def find_groups(measurements):
    """

    :param measurements:
    :return:
    """

    def build_graph(measurements):
        """

        :param measurements:
        :return:
        """
        graph = defaultdict(
            lambda: defaultdict(set))  # graph[i1][i2] = set(measurement idxs that connect these two iso_idxs)
        present = defaultdict(set)  # present[i] = set(measurements that touch this idx)

        for (m_idx, measurement) in measurements.iterrows():
            iso_idxs = measurement['Iso_Idxs'].split(',')
            for idx1, idx2 in product(iso_idxs, iso_idxs):
                if idx1 == idx2:
                    present[idx1].add(m_idx)
                else:
                    graph[idx1][idx2].add(m_idx)
                    graph[idx2][idx1].add(m_idx)
        return graph, present

    def graph_bfs(graph, present):
        """
        Explores graph that has an edge between idx1 and idx2 iff graph[idx1][idx2] exists


        :param graph:
            dict of dict of sets: graph[idx1][idx2] = set of all measurement indexes that include both isotopes
        :param present:
            dict of sets: present[idx] = set of all measurements that incldues this isotope
        :return:
            list of dicts: group['m_idxs'] is all measurement indexes in that group
                           group['iso_idxs'] is all isotopes in that group.
        """

        seen, groups = set(), []
        isotopes = set(present.keys())  # unique iso_idxs
        print(isotopes)
        for i in isotopes:
            if i in seen:
                continue

            group = {'m_idxs': set(), 'iso_idxs': set()}
            stack = [i]
            seen.add(i)
            while stack:
                curr = stack.pop()
                group['m_idxs'].update(present[curr])
                group['iso_idxs'].add(curr)

                for next_i in list(graph[curr].keys()):  # cant just do graph[curr] because we want to delete elements
                    group['m_idxs'].update(graph[curr][next_i])
                    del graph[curr][next_i]  # delete the edge of the graph
                    if next_i in seen:
                        continue
                    seen.add(next_i)  # mark as already visited
                    stack.append(next_i)

            group['m_idxs'] = list(group['m_idxs'])
            groups.append(group)

        return groups

    graph, present = build_graph(measurements)
    print(graph, present)
    return graph_bfs(graph, present)


def solve_group(group, measurements):
    print(group)
    params = lmfit.Parameters()

    for iso in group['iso_idxs']:
        a = int(iso[4:])
        params.add(f'R{iso}', 0.9 * a ** (1 / 3))

    for ridx, term in zip(group['m_idxs'], measurements.loc[group['m_idxs'], 'Term']):
        params.add(f'm{ridx}', expr=term)

    targets = list(measurements.loc[group['m_idxs'], ['Term', 'Value', 'Unc']].itertuples(index=False))
    targets += [
        (f"R{list(group['iso_idxs'])[0]}", 0, 100000)]  # add one absolute measurement to tether relative only groups

    def target(pars, tars):
        return [(pars.eval(term) - val) / unc for term, val, unc in tars]

    mini_res = lmfit.minimize(target, params, args=(targets,), scale_covar=False)

    return mini_res


def make_human_readable(measurements):
    measurements = measurements.copy()
    orig_columns = measurements.columns.tolist()

    num_isoidxs = measurements['Iso_Idxs'].str.split(',').apply(lambda i: len(i))
    max_idxs = max(num_isoidxs)
    measurements = measurements.reindex(columns=list(chain(*[[f'Z{i}', f'A{i}'] for i in range(1, max_idxs + 1)]))
                                                + orig_columns)

    for midx, row in measurements.iterrows():
        iso_idxs = row['Iso_Idxs'].split(',')
        for i, iso in enumerate(iso_idxs, 1):
            measurements.loc[midx, f'Z{i}'] = int(iso[:3])
            measurements.loc[midx, f'A{i}'] = int(iso[-3:])

    kept_columns = ['Value', 'Unc', 'Term', 'Connected_Measurements']
    measurements.drop(columns=[c for c in orig_columns if c not in kept_columns], inplace=True)

    absl = measurements.loc[num_isoidxs == 1, :]
    absl = absl.reindex(columns=['Z1', 'A1'] + [c for c in orig_columns if c != 'RefRemarks'])
    absl = absl.rename(columns={'Z1': 'Z', 'A1':'A'})

    rel = measurements.loc[num_isoidxs >= 2, :]

    return rel, absl


def save_human_readable(folder, measurements, group_info):
    rel_hr, abs_hr = make_human_readable(measurements)
    abs_hr.to_csv(f'./{folder}/absolute_radii_{folder}.csv', index=False)
    rel_hr.to_csv(f'./{folder}/relative_radii_{folder}.csv', index=False)
    group_info.to_csv(f'./{folder}/group_info_{folder}.csv', index=False)


if __name__ == '__main__':

    pd.set_option('display.width', None)
    meas_df = pd.read_csv('../measurements.csv')

    meas_absonly, group_info_absonly = compress_measurements(measurements=meas_df[meas_df['Table'].isin(['absl'])])
    save_human_readable('absonly', meas_absonly, group_info_absonly)

    meas_relonly, group_info_relonly = compress_measurements(
        measurements=meas_df[~meas_df['Table'].isin(['abs', 'noniso', 'ois'])])
    save_human_readable('relonly', meas_relonly, group_info_relonly)

