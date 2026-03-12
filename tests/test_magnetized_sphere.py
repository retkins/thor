"""Compute the magnetic field in a uniformly magnetized sphere

Sphere is 100mm diameter and centered at (0.0, 0.0, 1.0) meters in 3D space.
Sphere has mu_r = 1.2 and is a linear magnetic material (no B-H curve)
"""

from numpy.typing import NDArray
import pytest
import thor
import numpy as np
from numpy import float64
import matplotlib.pyplot as plt


@pytest.mark.parametrize(
    "method",
    ["dipole", "tet"],
)
def test_magnetized_sphere(method: str):
    # ---
    # Problem inputs
    # ---

    min_size: float = 10.0
    max_size: float = 20.0
    theta: float = 0.25
    alpha: float = 1.0  # relaxation factor for higher mu_r materials
    remesh: bool = True
    bext_mag = 5.0  # [T] External field magnitude
    mu_r = 1.2
    ctol: float = 1e-8  # convergence tolerance

    # ---
    # Load mesh and remesh if needed
    # ---

    # if remesh:
    #     thor.mesh.mesh_step("tests/data/sphere.stp", "tests/data/sphere_mesh.csv", min_size, max_size)
    # data = np.loadtxt("tests/data/sphere_mesh.csv", delimiter=',', skiprows=1)

    # centroids = data[:,0:3]
    # vol = data[:,3]
    nodes_flat, centroids, vol = thor.mesh.mesh_step_tets(
        "tests/data/sphere.stp", min_size, max_size
    )

    n = vol.shape[0]
    print("Thor - Magnetized sphere test")
    print(f"Problem size: {n} elements in mesh")

    # ---
    # Compute fields
    # ---

    # Apply a uniform H-field over the mesh
    h_ext = np.zeros((n, 3))
    h_ext[:, 2] = bext_mag / thor.MU0
    # h_ext[:,2] = np.linspace(0.0, bext_mag/thor.MU0, n)

    # Compute the magnetizing field acting on each element in the mesh
    # Total field is B = mu0 * (H + M)
    # M = (mu_r - 1) * H_ext
    mat = thor.LinearMaterial(mu_r)
    m_field = mat.chi(1) * h_ext
    moments = m_field * vol[:, np.newaxis]

    # Initial guess for h-field in the sphere using element dipole moments.
    if method == "dipole":
        h = thor.hfield_dipole(centroids, vol, moments, centroids, theta=theta)
    elif method == "tet":
        h = thor.hfield_dipole_tetrahedrons(
            nodes_flat,
            centroids,
            vol,
            m_field,
            centroids,
            theta=theta,
        )
    else:
        raise NotImplementedError("Unrecognized method")

    print(
        f"Initial guess for demagnetizing B-field in sphere (avg): {np.average(h, axis=0) * thor.MU0} T"
    )

    h_total = np.copy(h_ext)

    def fixed_point(h_total: NDArray[float64], n_iterations: int):
        """Compute fixed-point iteration on total field in the sphere

        Fixed point iteration: previous value of total field is used to calculate the demagntizing field, which updates the
        total field calculation.
        """

        for iteration in range(n_iterations):
            # Magnetization from current total field
            m_field = mat.chi(1) * h_total
            moments = m_field * vol[:, np.newaxis]  # m = M*V

            # Demagnetizing field from element dipole moments.
            if method == "dipole":
                h_demag = thor.hfield_dipole(
                    centroids, vol, moments, centroids, theta=theta
                )
            elif method == "tet":
                h_demag = thor.hfield_dipole_tetrahedrons(
                    nodes_flat,
                    centroids,
                    vol,
                    m_field,
                    centroids,
                    theta=theta,
                )
            else:
                raise NotImplementedError("Unrecognized method")

            # Update total field
            h_total_new = h_ext + h_demag

            # Check convergence
            change = np.max(np.abs(h_total_new - h_total))
            change_average = np.average(np.abs(h_total_new - h_total))
            max_idx = np.argmax(np.abs(change))
            print(f"Centroid: {centroids[max_idx]}")
            print(
                f"Distance from center: {np.linalg.norm(centroids[max_idx] - [0, 0, 1.0])}"
            )
            print(
                f"iter {iteration}: max change = {change:.6e}; avg change = {change_average:.6e}"
            )

            h_total = alpha * h_total_new + (1 - alpha) * h_total
            if change < ctol:
                print(f"Converged in {iteration + 1} iterations")
                print(f"\twith mu_r = {mu_r:.3f}")
                break

        return h_total

    h_total = fixed_point(h_total, 150)
    h_z_mean = np.mean(h_total[:, 2])
    h_z_analytical = (bext_mag / thor.MU0) / (1 + mat.chi(1) / 3)

    print(f"Mean H_z:       {h_z_mean:.0f} A/m")
    print(f"Analytical H_z: {h_z_analytical:.0f} A/m")
    print(f"Min B_z:        {np.min(h_total[:, 2]) * mu_r * thor.MU0:.3f} T")
    print(f"Max B_z:        {np.max(h_total[:, 2]) * mu_r * thor.MU0:.3f} T")
    print(f"Mean B_z:       {h_z_mean * mu_r * thor.MU0:.3f} T")
    print(f"Analytical B_z: {h_z_analytical * mu_r * thor.MU0:.3f} T")
    print(
        f"Error:          {abs(h_z_mean - h_z_analytical) / h_z_analytical * 100:.1f}%"
    )

    assert (np.abs(h_z_mean - h_z_analytical) / h_z_analytical) < 0.01

    b_z = h_total[:, 2] * mu_r * thor.MU0
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    s = ax.scatter(
        centroids[:, 0], centroids[:, 1], centroids[:, 2], c=b_z, cmap="viridis_r"
    )
    fig.colorbar(s)
    fig.savefig("tests/fig/mag_sphere.png")

    # Points within half an element size of the centerplane
    tol = max_size / 2 * 1e-3  # convert mm to meters
    mask = np.abs(centroids[:, 2] - 1.0) < tol

    fig, ax = plt.subplots()
    sc = ax.scatter(centroids[mask, 0], centroids[mask, 1], c=b_z[mask], cmap="viridis")
    plt.colorbar(sc, label="B_z [T]")
    ax.set_aspect("equal")
    ax.set_title("B_z on z-centerplane")
    plt.savefig("tests/fig/sphere_slice.png")
    plt.show()
