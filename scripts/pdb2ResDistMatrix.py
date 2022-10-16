# from a PDB file, calcuates the residue residue distance matrix
# for use in identifying residues that are closely contacting

# Oscar Charles 2022

import Bio.PDB
import numpy as np
import pandas as pd

pdb_code = ":)"
pdb_filename = "query/HCMV_UL54.pdb" 

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]

dist_matrix = calc_dist_matrix(model["A"], model["A"])
nloc = dist_matrix.shape[0] + 1


df = pd.DataFrame({ 'loc': range(1,nloc),
    'struc_mean_2_closest_residues': 0.0,
    'struc_mean_5_closest_residues': 0.0} )

for i in range(0 , dist_matrix.shape[0]) :
    print(i)
    # return the k closest residues
    k = 2
    dists = dist_matrix[i]
    idx = np.argpartition(dists, k+1)
    idx = np.delete(idx, 0) # remove self
    df.loc[i, 'struc_mean_2_closest_residues'] = np.mean(dists[idx[:k]])
    k = 5
    idx = np.argpartition(dists, k+1)
    idx = np.delete(idx, 0) # remove self
    df.loc[i, 'struc_mean_5_closest_residues'] = np.mean(dists[idx[:k]])

df.to_csv("/tmp/struc_mean_k_closest_residues.csv")