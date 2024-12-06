# Kuntal Ghosh
# May 2022

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.transform import Rotation

dr = 0.02
dang = 1.0
sigma = 10.0
r_fixed = 3.0
pi = np.pi

input_file = "../../cg.lammpstrj"
output_file = "adf.dat"

with open(input_file, "r") as f:
    lines = f.readlines()

# Extract number of atoms and box dimensions
natoms = int(lines[3].strip())
box_length = np.array([float(lines[5].split()[1]),
                       float(lines[6].split()[1]),
                       float(lines[7].split()[1])])

rho = natoms / np.prod(box_length)
nframes = len(lines) // (natoms + 9)  # 9 header lines per frame

# Allocate arrays for coordinates
xyz = np.zeros((nframes, natoms, 3))

# Parse coordinates
frame_idx = 0
line_idx = 0
while frame_idx < nframes:
    line_idx += 9  # Skip header lines
    for atom_idx in range(natoms):
        xyz[frame_idx, atom_idx, :] = np.array(
            list(map(float, lines[line_idx].split()[1:4]))
        )
        line_idx += 1
    frame_idx += 1

# Parameters for the angular distribution function (ADF)
angmin, angmax = 0.0, 180.0
nbin_ang = int((angmax - angmin) / dang) + 1
a = np.zeros(nbin_ang)

# Periodic boundary conditions
def apply_pbc(positions, box_length):
    """Apply periodic boundary conditions to ensure minimum image convention."""
    return positions - box_length * np.round(positions / box_length)

# Computes the angular distribution function
total = 0
for iframe in range(nframes):
    frame_coords = xyz[iframe]

    # Pairwise distances and vectors
    pairwise_vectors = frame_coords[:, None, :] - frame_coords[None, :, :]
    pairwise_vectors = apply_pbc(pairwise_vectors, box_length)
    distances = np.linalg.norm(pairwise_vectors, axis=-1)

    # Iterate over triplets of atoms
    for i in range(natoms):
        for j in range(natoms):
            if i != j:
                v1 = pairwise_vectors[i, j]
                d1 = distances[i, j]

                for k in range(j + 1, natoms):
                    if k != i:
                        v2 = pairwise_vectors[i, k]
                        d2 = distances[i, k]

                        # Check if distances are within the cutoff
                        if d1 < sigma and d2 < sigma:
                            # Calculate angle between vectors
                            cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                            ang = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

                            # Determine histogram bin
                            if angmin <= ang < angmax:
                                ibin_ang = int((ang - angmin) / dang)
                                a[ibin_ang] += 1
                                total += 1

# Normalize and write results
with open(output_file, "w") as f:
    for ibin_ang, count in enumerate(a):
        ang = angmin + ibin_ang * dang
        norm_value = count / (total * dang) if total > 0 else 0
        f.write(f"{ang * pi / 180.0:.6f} {norm_value:.6e}\n")

