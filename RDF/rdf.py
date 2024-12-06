# Kuntal Ghosh
# May 2022

import numpy as np
from scipy import integrate

# Constants
dr = 0.02
dang = 1.0
sigma = 3.7
r_fixed = 3.0
pi = np.pi

# Read the input file
with open('../traj.lammpstrj', 'r') as f:
    lines = f.readlines()

nlines = len(lines)
natoms = int(lines[3])
box_length = np.array([float(line.split()[1]) for line in lines[5:8]])

volume = np.prod(box_length)
rho = natoms / volume
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

dmin = 0.0
dmax = 0.5 * box_length[0]

print(f"dmin, dmax: {dmin}, {dmax}")

# Computing g(r)
nbin_r = int(0.5 * box_length[0] / dr) + 1
g = np.zeros(nbin_r)

for iframe in range(nframes):
    for i in range(natoms - 1):
        for j in range(i + 1, natoms):
            v = xyz[:, j, iframe] - xyz[:, i, iframe]
            v = v - box_length * np.round(v / box_length)
            d = np.linalg.norm(v)

            if 1.0 < d < (0.5 * box_length[0]):
                ibin_r = int((d - dmin) / dr)
                if 0 <= ibin_r < nbin_r:
                    g[ibin_r] += 2.0  # Count each pair twice

d_values = np.linspace(dmin, dmax, nbin_r)
shell_volume = 4 * np.pi * (d_values**2) * dr
ideal_count = rho * shell_volume * nframes * natoms

# Normalize g(r)
mask = ideal_count != 0
g[mask] /= ideal_count[mask]
g[~mask] = 0

# Save RDF data
np.savetxt('rdf.dat', np.column_stack((d_values, g)))

# Computing Nc (coordination number) using scipy's integrate function
d_values_nc = np.arange(dmin, sigma, dr)
integrand = 4.0 * pi * rho * (d_values_nc**2) * np.interp(d_values_nc, d_values, g)

# Filter out NaN and inf values
valid_indices = np.isfinite(integrand)
d_values_valid = d_values_nc[valid_indices]
integrand_valid = integrand[valid_indices]

nc = integrate.trapz(integrand_valid, d_values_valid)

print(f"Coordination number = {nc}")
