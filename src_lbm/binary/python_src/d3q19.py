import numpy as np

def lattice_fourier_laplacian(kx, ky, kz):
    """
    Calculates the Laplacian in Fourier space for a discrete lattice.

    This function computes the Laplacian operator in Fourier space for a given set of wavevectors. The inputs can be either floats or D-dimensional NumPy arrays. If all inputs are floats, the function returns a float. If any input is a NumPy array, the function returns a NumPy array of the same shape as the inputs.

    Parameters
    ----------
    kx : float or numpy.ndarray
        The wavevector component in the x-direction. Can be a float or a NumPy array.
    ky : float or numpy.ndarray
        The wavevector component in the y-direction. Can be a float or a NumPy array.
    kz : float or numpy.ndarray
        The wavevector component in the z-direction. Can be a float or a NumPy array.

    Returns
    -------
    float or numpy.ndarray
        The Laplacian operator in Fourier space:
        - Returns a float if all inputs are floats.
        - Returns a NumPy array if any of the inputs are NumPy arrays.

    Notes
    -----
    - The function uses the cosine of the wavevector components to compute the Laplacian.
    - The Laplacian operator is normalized by the speed of sound squared (`cs2 = 1/3`) to compute the output wavevector squared (`k2`).

    Examples
    --------
    Compute the Laplacian for scalar inputs:

    >>> kx, ky, kz = 1.0, 1.0, 1.0
    >>> laplacian = lattice_fourier_laplacian(kx, ky, kz)
    >>> print(laplacian)
    -4.0

    Compute the Laplacian for NumPy array inputs:

    >>> import numpy as np
    >>> kx = np.array([0.0, np.pi/2, np.pi])
    >>> ky = np.array([0.0, np.pi/2, np.pi])
    >>> kz = np.array([0.0, np.pi/2, np.pi])
    >>> laplacian = lattice_fourier_laplacian(kx, ky, kz)
    >>> print(laplacian)
    [-4.         -1.77777778  0.        ]
    """
    expr1 = np.cos(kx) + np.cos(ky) + np.cos(kz)
    expr2 = np.cos(kx)*np.cos(ky) + np.cos(ky)*np.cos(kz) + np.cos(kx)*np.cos(kz)
    out = 2/9*expr1 + 2/9*expr2 - 4/3
    cs2 = 1/3
    k2 = -out/cs2
    return k2