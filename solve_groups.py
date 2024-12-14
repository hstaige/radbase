import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

rel_df = pd.read_csv('relative_radii.csv')
abs_df = pd.read_csv('absolute_radii.csv')
iso_df = pd.read_csv('isotope_indices.csv')


def find_groups(iso, rel):
    graph = defaultdict(set)
    rdict = defaultdict(list)
    for r1, r2, ridx in zip(rel['iso_idx1'], rel['iso_idx2'], rel.index):
        graph[r1].add(r2)
        graph[r2].add(r1)

    isotopes = iso.index.tolist()
    seen = set()
    groups = []

    for i in isotopes:
        if i in seen:
            continue

        group = set()
        stack = [i]
        seen.add(i)
        while stack:
            curr = stack.pop()
            group.add(curr)

            for next_i in graph[curr]:
                if next_i in seen:
                    continue
                seen.add(next_i) #mark as already visited
                stack.append(next_i)
        if len(group) == 1:
            continue
        groups.append(group)

    return groups


groups = find_groups(iso_df, rel_df)
print([len(group) for group in groups])

f, ax = plt.subplots()

Z1, A1 = iso_df.loc[rel_df['iso_idx1']]['Z'].to_numpy(), iso_df.loc[rel_df['iso_idx1']]['A'].to_numpy()
Z2, A2 = iso_df.loc[rel_df['iso_idx2']]['Z'].to_numpy(), iso_df.loc[rel_df['iso_idx2']]['A'].to_numpy()

for z1, a1, z2, a2 in zip(Z1, A1, Z2, A2):
    c = 'k' if z1 == z2 else 'b'
    zorder = 0 if z1 == z2 else 10
    plt.plot([z1, z2], [a1-z1, a2-z2], f'-{c}s', zorder=zorder, markersize=3)

Z, A = iso_df.loc[abs_df['iso_idx']]['Z'].to_numpy(), iso_df.loc[abs_df['iso_idx']]['A'].to_numpy()
plt.plot(Z, A-Z, 'o', color='grey', zorder=-10, markersize=8)

plt.xlabel('Z')
plt.ylabel('N')
plt.show()