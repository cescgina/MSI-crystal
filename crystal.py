import htmd as ht
import numpy as np

filename = "1dxr.pdb"
f = open(filename, "r")
line = f.readline()
while(not line.startswith("EXPDTA")):
    line = f.readline()
if line.rfind("CRYSTAL") == -1:
    raise(ValueError("The input pdb does not correspond to a crystallographic "
                     "structure"))
while(not line.startswith("CRYST1")):
    line = f.readline()
(a, b, c, alpha, beta, gamma, group, g1, g2, g3, z) = line.split()[1:]
i = 0
orig_mat = np.zeros((3, 3))
T = np.zeros(3)
scale_mat = np.zeros((3, 3))
S = np.zeros(3)
line = f.readline()
while (i < 3):
    line_split = line.split()
    orig_mat[i, :] = line_split[1:4]
    T[i] = line_split[-1]
    i += 1
    line = f.readline()
i = 0
while (i < 3):
    line_split = line.split()
    scale_mat[i, :] = line_split[1:4]
    S[i] = line_split[-1]
    i += 1
    line = f.readline()
f.close()
mol = ht.Molecule(filename)
# http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
# http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction
