from collections import defaultdict
from itertools import chain, product

import lmfit
import numpy as np
import pandas as pd
from uncertainties import wrap as uwrap


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
    return graph_bfs(graph, present)


def solve_group(group, measurements):

    def target(pars, tars):
        # note, the two zeros are to force scalar minimize to calc covariances ...
        return [0] * 1 + [(pars.eval(term) - val) / unc for term, val, unc in tars]

    def iter_cb(params, iter, resids, *args, **kwargs):
        nfree = len([p for p in params if params[p].vary])
        if nfree < 50 or (iter % 100) != 0:
            return
        print(f'On iter: {iter} with nfree = {nfree}')
        print(f'Current redchi is {(np.array(resids) ** 2).sum() / (len(resids) - nfree)}')
        print(sorted(list(zip(resids, targets[2:])), reverse=True)[:10])

    def calc_measurement_uncs(minres, measure, grp):
        # addresses the issue of calculating uncertainties being very expensive when there are a large number of parameters
        # by fixing parameters that do not affect that quantity. E.G. if m1 = 'R001002-R001001' then we don't need to vary
        # R060120 when calculating the uncertainty of m1.

        # Essentially copies lmfit.parameters.create_uvars

        wrap_ueval = uwrap(lmfit.parameter.asteval_with_uncertainties)
        pars = minres.params
        vnames = list(pars.keys())
        for ridx, (term, iso_idx) in zip(grp['m_idxs'], measure.loc[grp['m_idxs'], ['Term', 'Iso_Idxs']].itertuples(index=False)):
            mname = f'm{ridx}'
            pars.add(mname, expr=term)

            relevant = ['R' + idx for idx in iso_idx.split(',')]
            corr_vars = [minres.uvars[p] if p in relevant else pars[p].value for p in vnames]
            uval = wrap_ueval(*corr_vars, obj=pars[mname], pars=pars, names=vnames)
            pars[mname].stderr = uval.std_dev
            minres.uvars[mname] = uval

        minres.residual = minres.residual[1:]  # drop extra 0 for scalar minimizers
        minres._calculate_statistics()

        return minres

    params = lmfit.Parameters()

    for iso in group['iso_idxs']:
        a = int(iso[3:])
        params.add(f'R{iso}', (0.9071 + 1.1025 / (a ** (2 / 3)) + -0.548 / (a ** (4 / 3))) * a ** (1 / 3))

    targets = list(measurements.loc[group['m_idxs'], ['Term', 'Value', 'Unc']].itertuples(index=False))

    method = 'leastsq'
    mini = lmfit.Minimizer(target, params, fcn_args=(targets,), scale_covar=False, calc_covar=True, iter_cb=iter_cb,
                           max_nfev=100000)

    mini_res = mini.minimize(method=method)

    if not mini_res.errorbars and len(group['iso_idxs']) > 1:
        params = mini_res.params
        params[f'R{list(group['iso_idxs'])[0]}'].vary = False
        mini_res = mini.minimize(params=params, method=method)
        if not mini_res.errorbars:
            raise AttributeError('No errorbars for group', group)

    mini_res = calc_measurement_uncs(mini_res, measurements, group)

    return mini_res


def make_human_readable(measurements):
    measurements = measurements.copy()
    orig_columns = measurements.columns.tolist()

    num_isoidxs = measurements['Iso_Idxs'].str.split(',').apply(lambda i: len(i))
    max_idxs = max(num_isoidxs)
    measurements = measurements.reindex(columns=list(chain(*[[f'Z{i}', f'A{i}'] for i in range(1, max_idxs + 1)])) + orig_columns)

    for midx, row in measurements.iterrows():
        iso_idxs = row['Iso_Idxs'].split(',')
        for i, iso in enumerate(iso_idxs, 1):
            measurements.loc[midx, f'Z{i}'] = int(iso[:3])
            measurements.loc[midx, f'A{i}'] = int(iso[-3:])

    kept_columns = ['Value', 'Unc', 'Term', 'Connected_Measurements']
    measurements.drop(columns=[c for c in orig_columns if c not in kept_columns], inplace=True)

    absl = measurements.loc[num_isoidxs == 1, :]
    absl = absl.reindex(columns=['Z1', 'A1'] + [c for c in orig_columns if c != 'RefRemarks'])
    absl = absl.rename(columns={'Z1': 'Z', 'A1': 'A'})

    rel = measurements.loc[num_isoidxs >= 2, :]

    return rel, absl


def save_human_readable(folder, measurements, group_info, suffix=''):
    rel_hr, abs_hr = make_human_readable(measurements)
    abs_hr.to_csv(f'./{folder}/absolute_radii_{suffix}.csv', index=False)
    rel_hr.to_csv(f'./{folder}/relative_radii_{suffix}.csv', index=False)
    group_info.to_csv(f'./{folder}/group_info_{suffix}.csv', index=False)


if __name__ == '__main__':
    pd.set_option('display.width', None)
    meas_df = pd.read_csv('../inputs/measurements.csv')

    print('---- Abs only ----')
    meas_absonly, group_info_absonly = compress_measurements(measurements=meas_df[meas_df['Table'].isin(['absl'])])
    save_human_readable('../outputs/absonly', meas_absonly, group_info_absonly, 'absonly')

    print('---- rel noniso only ----')
    meas_relonly, group_info_relonly = compress_measurements(
        measurements=meas_df[~meas_df['Table'].isin(['abs', 'noniso', 'ois', 'kalpha'])])
    save_human_readable('../outputs/relonly', meas_relonly, group_info_relonly, 'relonly')

    print('---- All ----')
    meas_relonly, group_info_relonly = compress_measurements(
        measurements=meas_df[~meas_df['Table'].isin(['kalpha'])])
    save_human_readable('../outputs/all_nokalpha', meas_relonly, group_info_relonly, 'all_nokalpha')
