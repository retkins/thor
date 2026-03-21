"""Compute force on magnetic material (uniform sphere)"""

import numpy as np
import oersted

min_size: float = 15
max_size: float = 25
b_ext: float = 5.0
material = oersted.LinearMaterial(1.5)

nodes, centroids, vol = oersted.mesh.mesh_step_tets("tests/data/sphere.stp", min_size, max_size)

# create a uniform gradient across the mesh in z-direction
z_min = np.min(centroids[:, 2])
z_max = np.max(centroids[:, 2])
h0 = b_ext / oersted.MU0
dhdz = h0 / (z_max - z_min)


def h_ext(targets):
    h = np.zeros(targets.shape)
    h[:, 2] = h0 / 2.0 + h0 * (targets[:, 2] - z_min) / (z_max - z_min)
    return h


forces = oersted.mag_force(centroids, vol, material, h_ext)

fmag = np.linalg.norm(np.sum(forces, axis=0))
f_expected = oersted.MU0 * material.chi(1) * np.sum(vol) * h0 * dhdz
err = (f_expected - fmag) / fmag

print("oersted - magnetized material test\n---\n")
print(f"total force:    {fmag:.0f} N")
print(f"expected force: {f_expected:.0f} N")
print(f"error: {100 * err:.3f} %")


def test_mag_force():
    assert np.abs(err) < 0.01
