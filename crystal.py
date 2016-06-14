import htmd as ht
import numpy as np
import copy
# import matplotlib.pyplot as plt

# from mpl_toolkits.mplot3d import Axes3D


class MoleculeCopy():

    def __init__(self, mol):
        self.mol = copy.deepcopy(mol)

    def is_outside(self, size):
        center = np.mean(self.mol.get('coords'), axis=0)
        cellloc = center / size
        cellloc = np.floor(cellloc)
        self.cellloc = cellloc
        self.center = center
        self.box_size = size
        return np.sum(np.abs(cellloc)) > 0

    def move_inside(self, axes):
        neworigin = np.multiply(axes.transpose(), self.cellloc).transpose()
        neworigin = np.sum(neworigin, axis=0)
        newcenter = self.center-neworigin
        self.mol.center(loc=newcenter)


class UnitCell():

    def __init__(self, pdbfile):
        self.mol = ht.Molecule(pdbfile)
        self.readPDB(pdbfile)
        caux = (np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma)) / \
            np.sin(self.gamma)
        self.axes = np.array([[self.a, 0, 0], [self.b*np.cos(self.gamma),
                                               self.b*np.sin(self.gamma), 0],
                              [self.c*np.cos(self.beta), self.c * caux,
                               self.c * np.sqrt(1-np.cos(self.beta)**2-caux**2)
                               ]])
        self.molunit = ht.Molecule()
        self.viewer = ht.vmdviewer.VMD()
        self.draw_cell()
        self.num_copies = len(self.trans_v)
        self.draw_copies()

    def readPDB(self, pdbfile):

        f = open(filename, "r")
        line = f.readline()
        while(not line.startswith("EXPDTA")):
            line = f.readline()
        if (line.rfind("CRYSTAL") == -1) and (line.rfind("DIFFRACTION") == -1):
            raise(ValueError("The input pdb does not correspond to a "
                             "crystallographic structure"))
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
        self.rot_mat = rot_mat
        self.trans_v = trans_v
        while(not line.startswith("CRYST1")):
            line = f.readline()
        (a, b, c, alpha, beta, gamma) = line.split()[1:7]
        self.group = line.split()[7:-2]
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
        self.orig_mat = orig_mat
        self.scale_mat = scale_mat
        self.T = T
        self.S = S
        self.alpha = np.deg2rad(float(alpha))
        self.beta = np.deg2rad(float(beta))
        self.gamma = np.deg2rad(float(gamma))
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)

    def draw_cell(self):

        baseorigin = np.array([0, 0, 0])
        base_lo_l = self.axes[0]
        base_up_r = self.axes[1]
        base_lo_r = self.axes[0] + self.axes[1]
        top_origin = self.axes[2]
        top_up_r = self.axes[1] + self.axes[2]
        top_lo_l = self.axes[0] + self.axes[2]
        top_lo_r = self.axes.sum(axis=0)
        points = np.vstack((baseorigin, base_lo_l, base_up_r, base_lo_r))
        points = np.vstack((points, top_origin, top_lo_l, top_up_r, top_lo_r))
        size = np.array([base_lo_l[0], base_up_r[1], top_origin[2]])
        if self.group[0] == 'H':
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
            lines_draw = np.array([[1, 3], [1, 5], [2, 3], [2, 6], [3, 7],
                                   [5, 7], [6, 7]])
            lines_draw = np.vstack((lines_draw, lines_draw+8, lines_draw+16))
        else:
            lines_draw = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3],
                          [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]]
        for row in lines_draw:
            self.viewer.send('draw line {{{0} {1} {2}}} {{{3} {4} {5}}}'.format(
                *(tuple(points[row[0]]) + tuple(points[row[1]]))))
        self.size = size

    def draw_copies(self):
        for i in range(self.num_copies):
            molecule = MoleculeCopy(self.mol)
            molecule.mol.rotateBy(self.rot_mat[i])
            molecule.mol.moveBy(self.trans_v[i])
            if molecule.is_outside(self.size):
                molecule.move_inside(self.axes)
            self.molunit.append(molecule.mol)

    def view_unit_cell(self):
        self.molunit.view(style='NewCartoon', viewerhandle=self.viewer)

if __name__ == "__main__":
    filename = "3v03.pdb"
    unitcell = UnitCell(filename)
    unitcell.view_unit_cell()
