import htmd as ht
import numpy as np

filename = "1dxr.pdb"
f = open(filename, "r")
line = f.readline()
while(not line.startswith("EXPDTA")):
    line = f.readline()
if (line.rfind("CRYSTAL") == -1) and (line.rfind("DIFFRACTION") == -1):
    raise(ValueError("The input pdb does not correspond to a crystallographic "
                     "structure"))
# find start crystal structures
while(line.rfind("REMARK 290       1555") == -1):
    line = f.readline()
# count structures
nstr = 0
while(len(line.split())>2):
	nstr += 1
	print(len(line.split()))
	line = f.readline()
rot_mat = np.zeros((nstr, 3, 3)) # generate empty coord matrix
trans_v = np.zeros((nstr, 3))

while(line.rfind("SMTRY1") == -1):
	line = f.readline()
	
# fill rot_mat and trans_v
n = 0
dim = 0
while(len(line.split())>2):	
	values = line.split()
	rot_mat[n][dim] = values[4:7]
	print(rot_mat[n][dim])
	trans_v[n][dim] = values[7]
	print(trans_v[n])
	line = f.readline()
	if dim == 2:
		n += 1
		dim = 0
	else:
		dim += 1
print(rot_mat)
print(trans_v)

mol = ht.Molecule(filename)

# http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
# http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction
