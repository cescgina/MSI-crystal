import htmd as ht
import numpy as np
# import matplotlib.pyplot as plt
# import copy
# from mpl_toolkits.mplot3d import Axes3D


def pack_mol(mol, max_vec=None):
    old_shape = mol.coords.shape
    traj_flat = mol.coords.transpose(2, 0, 1).reshape(-1, mol.coords.shape[1])
    traj_flat = np.mod(traj_flat, mol.box.transpose())
    mol.coords = traj_flat.reshape(old_shape)
    return mol


def draw_cell(axis, group, viewer):

    baseorigin = np.array([0, 0, 0])
    base_lo_l = axis[0]
    base_up_r = axis[1]
    base_lo_r = axis[0] + axis[1]
    top_origin = axis[2]
    top_up_r = axis[1] + axis[2]
    top_lo_l = axis[0] + axis[2]
    top_lo_r = axis.sum(axis=0)
    points = np.vstack((baseorigin, base_lo_l, base_up_r, base_lo_r))
    points = np.vstack((points, top_origin, top_lo_l, top_up_r, top_lo_r))
    if group[0] == 'H':
        # Lattice is hexagonal, can be constructed by applying two 120º
        # rotations over the z axis
        alpha = np.deg2rad(120)
        rot_mat = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])
        rotpoints1 = np.dot(rot_mat, points.transpose())
        rotpoints2 = np.dot(rot_mat, rotpoints1)
        points = np.vstack((points, rotpoints1.transpose(),
                           rotpoints2.transpose()))
        lines_draw = [[1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [6, 7]]
    else:
        lines_draw = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6],
                      [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]]
    for row in lines_draw:
        viewer.send('draw line {{{0} {1} {2}}} {{{3} {4} {5}}}'.format(
            *(tuple(points[row[0]]) + tuple(points[row[1]]))))

filename = "3v03.pdb"
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
while(len(line.split()) > 2):
    nstr += 1
    line = f.readline()
rot_mat = np.zeros((nstr, 3, 3))  # generate empty coord matrix
trans_v = np.zeros((nstr, 3))

while(line.rfind("SMTRY1") == -1):
    line = f.readline()

# fill rot_mat and trans_v
n = 0
dim = 0
while(len(line.split()) > 2):
    values = line.split()
    rot_mat[n][dim] = values[4:7]
    trans_v[n][dim] = values[7]
    line = f.readline()
    if dim == 2:
        n += 1
        dim = 0
    else:
        dim += 1

while(not line.startswith("CRYST1")):
    line = f.readline()
(a, b, c, alpha, beta, gamma) = line.split()[1:7]
group = line.split()[7:-2]
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
alpha = np.deg2rad(float(alpha))
beta = np.deg2rad(float(beta))
gamma = np.deg2rad(float(gamma))
a = float(a)
b = float(b)
c = float(c)
mol = ht.Molecule(filename)
caux = (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
axes = np.array([[a, 0, 0], [b*np.cos(gamma), b*np.sin(gamma), 0],
                [c*np.cos(beta), c * caux,
                 c * np.sqrt(1-np.cos(beta)**2-caux**2)]])
c2r = axes.transpose()
r2c = np.linalg.inv(c2r)
asum = np.add(axes[0], axes[1])
np.add(asum, axes[2], asum)
origin = np.array([0, 0, 0])
center = np.add(origin, 0.5*asum)
mol2 = ht.Molecule()
viewer = ht.vmdviewer.VMD()
for i in range(len(trans_v)):
    molecule = mol.copy()
    molecule.rotateBy(rot_mat[i])
    molecule.moveBy(trans_v[i])
    molecule.view(style='New Cartoon', viewerhandle=viewer)
    mol2.append(molecule)
# mol2.view()
draw_cell(axes, group, viewer)
# http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
# http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction
