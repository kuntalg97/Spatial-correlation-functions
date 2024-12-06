# Kuntal Ghosh
# Code for computing adngular distribution functions (ADFs)
# May 2022

import numpy as np

# Constants
dr = 0.02
dang = 1.0
sigma = 10.0
r_fixed = 3.0
pi = np.pi

def angle(v1, v2):
    cos_theta = np.dot(v1, v2)
    d1 = np.linalg.norm(v1)
    d2 = np.linalg.norm(v2)
    cos_theta = cos_theta / (d1 * d2)
    a = np.arccos(cos_theta)
    return np.degrees(a)

# Read the input file
with open('../traj.lammpstrj', 'r') as f:
    lines = f.readlines()

nlines = len(lines)
natoms = int(lines[3])
box_length = np.array([float(line.split()[1]) for line in lines[5:8]])
print(box_length)

rho = natoms / np.prod(box_length)
nframes = nlines // (natoms + 9)

# Initialize arrays
xyz = np.zeros((3, natoms, nframes))
force_xyz = np.zeros((3, natoms, nframes))

# Store coordinates
for i in range(nframes):
    start = i * (natoms + 9) + 9
    for k in range(natoms):
        line = lines[start + k].split()
        xyz[:, k, i] = [float(x) for x in line[1:4]]
        force_xyz[:, k, i] = [float(x) for x in line[4:7]]

angmin, angmax = 0.0, 180.0
dmin, dmax = 0.0, 0.5 * box_length[0]

print(f"angmin, angmax, dmin, dmax: {angmin}, {angmax}, {dmin}, {dmax}")

nbin_ang = int((angmax - angmin) / dang) + 1
a = np.zeros(nbin_ang)

total = 0
for iframe in range(nframes):
    for i in range(natoms):
        for j in range(natoms):
            if i != j:
                v1 = xyz[:, j, iframe] - xyz[:, i, iframe]
                v1 -= box_length * np.round(v1 / box_length)
                d1 = np.dot(v1, v1)

                for k in range(j + 1, natoms):
                    if k != i:
                        v2 = xyz[:, k, iframe] - xyz[:, i, iframe]
                        v2 -= box_length * np.round(v2 / box_length)
                        d2 = np.dot(v2, v2)

                        if d1 < sigma and d2 < sigma:
                            ang = angle(v1, v2)

                            if 0.0 < ang < 180.0:
                                ibin_ang = int((ang - angmin) / dang)

                                if 0 <= ibin_ang < nbin_ang:
                                    a[ibin_ang] += 1.0
                                    total += 1

# Save ADF data
ang_values = np.arange(angmin, angmax + dang, dang)
adf = a / (total * dang)
np.savetxt('adf.dat', np.column_stack((np.radians(ang_values), adf)))
