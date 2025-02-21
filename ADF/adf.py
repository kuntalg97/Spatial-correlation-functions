# Kuntal Ghosh
# Code for computing three-body correlation functions
# Modify the types for J-I-K type triplets in this code

import numpy as np
from numba import jit
import matplotlib.pyplot as plt

@jit(nopython=True)
def angle(v1, v2):
    cos_theta = np.dot(v1, v2)
    d1 = np.linalg.norm(v1)
    d2 = np.linalg.norm(v2)
    cos_theta /= (d1 * d2)
    return np.arccos(cos_theta) * 180 / np.pi

@jit(nopython=True)
def compute_adf(xyz, type_id, box_length, nframes, natoms, dr, dang, sigma, angmin, angmax, target_i, target_j, target_k):
    nbin_ang = int((angmax - angmin) / dang) + 1
    a = np.zeros(nbin_ang)
    total = 0

    for iframe in range(nframes):
        for i in range(natoms):
            if type_id[i] != target_i:
                continue
            for j in range(natoms):
                if i == j or type_id[j] != target_j:
                    continue
                for k in range(j + 1, natoms):
                    if k == i or type_id[k] != target_k:
                        continue

                    v1 = xyz[j, :, iframe] - xyz[i, :, iframe]
                    v1 -= box_length * np.round(v1 / box_length)
                    d1 = np.linalg.norm(v1)

                    v2 = xyz[k, :, iframe] - xyz[i, :, iframe]
                    v2 -= box_length * np.round(v2 / box_length)
                    d2 = np.linalg.norm(v2)

                    if d1 < sigma and d2 < sigma:
                        ang = angle(v1, v2)
                        ibin_ang = int((ang - angmin) / dang)
                        if 0 <= ibin_ang < nbin_ang:
                            a[ibin_ang] += 1
                            total += 1

    return a, total

def main():
    # Parameters
    dr = 0.02
    dang = 0.5
    sigma = 6.5
    angmin, angmax = 0.0, 180.0

    # Target atom types for J-I-K triplet
    target_j, target_i, target_k = map(int, input("Enter the atom types for J, I (central atom), K: ").split())

    # Read LAMMPS trajectory file
    with open('../cg_np_wat.lammpstrj', 'r') as f:
        lines = f.readlines()

    natoms = int(lines[3])
    box_length = np.array([float(line.split()[1]) - float(line.split()[0]) for line in lines[5:8]])
    nframes = len(lines) // (natoms + 9)

    xyz = np.zeros((natoms, 3, nframes))
    force_xyz = np.zeros((natoms, 3, nframes))
    type_id = np.zeros(natoms)

    for i in range(nframes):
        frame_start = i * (natoms + 9) + 9
        frame_end = frame_start + natoms
        frame_data = np.array([line.split() for line in lines[frame_start:frame_end]], dtype=float)
        type_id = frame_data[:, 1]
        xyz[:, :, i] = frame_data[:, 2:5]
        force_xyz[:, :, i] = frame_data[:, 5:8]

    # Compute ADF
    a, total = compute_adf(xyz, type_id, box_length, nframes, natoms, dr, dang, sigma, angmin, angmax, target_i, target_j, target_k)

    # Write results
    nbin_ang = int((angmax - angmin) / dang) + 1
    ang_range = np.linspace(angmin, angmax, nbin_ang)
    adf = a / (total * dang)

    out_file = f'adf_{target_j}{target_i}{target_k}.dat'
    np.savetxt(out_file, np.column_stack((ang_range, ang_range * np.pi / 180, adf)))

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(ang_range, adf)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Probability')
    plt.title(f'Angular Distribution Function (J-I-K: {target_j}-{target_i}-{target_k})')
    plt.savefig('adf_plot.png')
    plt.show()

if __name__ == "__main__":
    main()

