
from collections import Counter, defaultdict
from ase.neighborlist import neighbor_list
import numpy as np

def get_neighbor_count(structure,cutoffs):
        cutoffs.sort(reverse=True)
        nth_shell_counts = []
        for cutoff in cutoffs:
            i, j = neighbor_list("ij", structure.atoms, cutoff)
            #self._total_bond.append(len(i))
            data = defaultdict(list)
            for p, q in zip(i, j):
                data[structure.atom_dict[p]].append(structure.atom_dict[q])
            nth_shell_counts.append({key: Counter(data[key]) for key in data.keys()})

        for idx in range(len(nth_shell_counts) - 1):
            for key in nth_shell_counts[idx].keys():
                nth_shell_counts[idx][key].subtract(nth_shell_counts[idx + 1][key])
        return nth_shell_counts

def get_neighbor_count_matrix(ndatas):
    ndata_matrices=[]
    for ndata in ndatas:
        elements=ndata.keys()
        datas=[]
        for element_A in elements:
            data=[]
            for element_B in elements:
                data.append(ndata[element_A][element_B])
            datas.append(data)
        ndata_matrices.append(datas)

    return np.array(ndata_matrices)