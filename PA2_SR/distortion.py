"""
Bernstein Polynomial Distortion Correction for PA2

This module implements distortion correction using Bernstein polynomials.

Sidharth Raghavan
"""

import numpy as np
from typing import Tuple
from scipy.special import comb


def bernstein_basis(u: np.ndarray, i: int, N: int) -> np.ndarray:
    """
    Calculate Bernstein basis polynomial B_{N,i}(u).
    B_{N,i}(u) = C(N,i) * u^i * (1-u)^{N-i}; C(N,i) is the binomial coefficient
    """
    #handle edge cases
    u = np.clip(u, 0, 1) 
    
    #binomial coefficient
    c_ni = comb(N, i, exact=True)
    
    #bernstein polynomial
    return c_ni * np.power(u, i) * np.power(1 - u, N - i)


def get_index(i: int, j: int, k: int, N: int) -> int:
    """
    Calculate index for Bernstein polynomial coefficient.
    """
    return int(i + j * (N + 1) + k * (N + 1)**2)


def bernstein_polynomial(u: np.ndarray, N: int = 5) -> np.ndarray:
    """
    compute Bernstein polynomial basis matrix F.
    
    for each point (x, y, z), compute
    F_ijk(x, y, z) = B_{N,i}(x) * B_{N,j}(y) * B_{N,k}(z)

    return the basis matrix F with shape N_points, (N+1)^3
    """
    N_points = u.shape[0]
    n_coeffs = (N + 1)**3
    
    F = np.zeros((N_points, n_coeffs))
    
    for point_idx in range(N_points):
        ux, uy, uz = u[point_idx, 0], u[point_idx, 1], u[point_idx, 2]
        
        #get all comb of i, j, k
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    # get the individual Bernstein basis polynomials
                    Bx = bernstein_basis(ux, i, N)
                    By = bernstein_basis(uy, j, N)
                    Bz = bernstein_basis(uz, k, N)
                    
                    #combined basis
                    index = get_index(i, j, k, N)
                    F[point_idx, index] = Bx * By * Bz
    
    return F


def scale_to_box(q: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """
    Scale data to [0, 1] box while preserving proportions
    """
    qmin = np.min(q)
    qmax = np.max(q)
    
    # Avoid division by zero
    if abs(qmax - qmin) < 1e-10:
        return q, qmin, qmax
    
    v = (q - qmin) / (qmax - qmin)
    return v, qmin, qmax


def unscale_from_box(v: np.ndarray, q: np.ndarray = None, qmin: float = None, qmax: float = None) -> np.ndarray:
    """
    Rescale data from [0, 1] box back to original range
    """
    if q is not None:
        qmin = np.min(q)
        qmax = np.max(q)
    
    if qmin is None or qmax is None:
        raise ValueError("need either q or both qmin and qmax")
    
    return v * (qmax - qmin) + qmin


def distortion_coeff(true_points: np.ndarray, distorted_points: np.ndarray, N: int = 5) -> np.ndarray:
    """
    Find distortion coefficients using Bernstein polynomials.

    Input: ground truth points, distorted 3D points, polynomial degree
    output- bernstein polynomial coefficients
    """
    #included for case checking
    if len(true_points.shape) == 3:
        true_points = np.reshape(np.transpose(true_points, (0, 2, 1)), (-1, 3))
    if len(distorted_points.shape) == 3:
        distorted_points = np.reshape(np.transpose(distorted_points, (0, 2, 1)), (-1, 3))
    
    #scale data
    true_scaled, _, _ = scale_to_box(true_points)
    distorted_scaled, _, _ = scale_to_box(distorted_points)
    
    #compute matrix
    F = bernstein_polynomial(distorted_scaled, N)
    
    #get coefficients
    from scipy.linalg import lstsq
    coeff, _, _, _ = lstsq(F, true_scaled)
    
    return coeff


def distortion_correction(distorted_values: np.ndarray, coeff: np.ndarray, N: int = 5) -> np.ndarray:
    """
    correct distortion using bernstein polynomial coefficients

    Outputs the corrected points with same shape as distorted_values
    """
    og_shape = distorted_values.shape
    is_3d = len(og_shape) == 3
    
    if is_3d:
        distorted_values_2d = np.reshape(np.transpose(distorted_values, (0, 2, 1)), (-1, 3))
    else:
        distorted_values_2d = distorted_values
    
    #scale
    val_scaled, _, _ = scale_to_box(distorted_values_2d)
    
    #get basis matrix
    F = bernstein_polynomial(val_scaled, N)
    
    #correction
    corrected_scaled = F @ coeff
    
    #rescale
    corrected_flat = unscale_from_box(corrected_scaled, q=distorted_values_2d)
    
    if is_3d:
        corrected = np.transpose(np.reshape(corrected_flat, (og_shape[0], og_shape[2], og_shape[1])), (0, 2, 1))
    else:
        corrected = corrected_flat
    
    return corrected


def main():
    print("Testing Bernstein polynomial distortion correction...")
    
    #test points 
    true_points = np.array([[0.0, 0.0, 0.0], [100.0, 100.0, 100.0], [50.0, 50.0, 50.0]])
    
    #manually add distortion
    distorted_points = true_points + np.random.normal(0, 0.1, true_points.shape)
    
    print(f"True points:\n{true_points}")
    print(f"Distorted points:\n{distorted_points}")
    
    coeff = distortion_coeff(true_points, distorted_points)
    print(f"\nCoefficients computed: shape {coeff.shape}")
    
    corrected = distortion_correction(distorted_points, coeff)
    print(f"Corrected points:\n{corrected}")
    
    error = np.linalg.norm(corrected - true_points)
    print(f"\nReconstruction error: {error:.6f}")


if __name__ == "__main__":
    main()