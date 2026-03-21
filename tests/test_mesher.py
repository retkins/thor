"""Confirm that the mesher is working properly

Sphere is of radius r=50mm
"""

import oersted
import numpy as np
from time import perf_counter

min_size: float = 5.0
max_size: float = 5.0

start = perf_counter()
nodes, centroids, volumes = oersted.mesh.mesh_step_tets("tests/data/sphere.stp", min_size, max_size)
end = perf_counter()
n = len(volumes)


r: float = 50.0
vol_expected = (4.0 / 3.0) * np.pi * 50**3
vol_mesh: float = np.sum(volumes) * 1e9  # convert to mm3
vol_err: float = (vol_mesh - vol_expected) / vol_expected
print(f"Meshed {n} tet elements in {end - start:.3f} seconds")
print(f"Total volume: \n\tmesh:   {vol_mesh:.0f} mm3 \n\tsphere: {vol_expected:.0f} mm3")
print(f"\tdiff = {vol_err * 100:.2f} %")


def test_mesher():
    assert np.abs(vol_err) < 0.01


if __name__ == "__main__":
    test_mesher()
