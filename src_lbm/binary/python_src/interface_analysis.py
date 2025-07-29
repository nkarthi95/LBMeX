import numpy as np
from skimage.measure import find_contours
from swift_orlandini import swift_theoretical_C0, swift_theoretical_xi
from scipy.optimize import curve_fit

def bin_height_spectrum_1D(k1, S1, nbins): 
    """
    Bin a 1D Fourier spectrum by averaging intensities over specified frequency bins.

    Parameters
    ----------
    k1 : np.ndarray
        A 1D numpy array of frequencies obtained from a Fourier transform.
    S1 : np.ndarray
        A 1D numpy array representing the intensities of the data in the Fourier transform.
    nbins : int
        The number of bins to use for binning the data.

    Returns
    -------
    xbin : np.ndarray
        A 1D numpy array of binned frequency values.
    ybin : np.ndarray
        A 1D numpy array of binned intensity values.
    
    Notes
    -----
    - The function selects the first half of the spectrum (excluding the DC component).
    - Binning is performed using `np.histogram`, averaging intensity values within each bin.
    - The bin centers are computed as the midpoint of each bin's edges.
    """

    L = S1.size
    slc = slice(1, L//2-1)
    xraw = k1[slc].copy()
    yraw = S1[slc].copy()

    kmin = 2*np.pi/nbins
    bins = np.arange(nbins//2+1)*kmin # kmax+1 for bin_edges: len(bins)=len(hist)+1
    shells = np.histogram(xraw, bins, weights=yraw)[0]
    counts = np.histogram(xraw, bins)[0]

    xbin = (bins[:-1]+bins[1:])/2
    ybin = shells/counts

    return xbin, ybin

def ih_direct(profile, level = 0):
    """
    Calculates the interface height in the x-direction for a system with a flat interface across the yz plane.

    This function determines the x-coordinate of the interface at each point in the yz plane by identifying where the isocontour of the interface is zero. It uses the `skimage.measure.find_contours` function to perform this calculation.

    Parameters
    ----------
    profile : numpy.ndarray
        A 3D NumPy array representing the order parameter of the system. The array has dimensions `(nx, ny, nz)`, where `nx` is the size in the x-direction and `ny` and `nz` are the sizes in the y- and z-directions, respectively.

    Returns
    -------
    numpy.ndarray
        A 2D NumPy array of shape `(ny, nz)` containing the x-heights of the interface for each point in the yz plane.

    Notes
    -----
    - The interface is assumed to be flat across the yz plane.
    - The `skimage.measure.find_contours` function is used to locate the zero isocontour of the interface along the x-direction for each slice in the yz plane.

    Examples
    --------
    Calculate the interface height for a given profile:

    >>> import numpy as np
    >>> from skimage import measure
    >>> profile = np.random.random((50, 30, 30)) - 0.5  # Example 3D order parameter
    >>> interface_heights = ih_direct(profile)
    >>> print(interface_heights.shape)
    (30, 30)
    """
    nx, ny, nz = profile.shape
    index_levels = np.arange(0, nz, 1)
    out = np.zeros((ny, nz))
    for y in range(ny):
        levels = find_contours(profile[:, y, :], level = level)[0]
        idxs = np.isin(levels[:, 1], index_levels)
        out[y] = levels[idxs, 0]
    return out

def ih_profile_fit(profile, chi, T, kappa):
    nx, ny, nz = profile.shape
    out = np.zeros((ny, nz))

    xidxs = np.arange(0, nx, 1)
    slc = slice(nx//4, 3*nx//4)
    xfit = xidxs[slc]
    
    ce = swift_theoretical_C0(chi/T)
    phi0 = 2*ce - 1
    zeta = swift_theoretical_xi(ce, chi/T, kappa)

    fit_func = lambda x, hint, p0, w: p0*np.tanh((x - hint)/(w))

    for y in range(ny):
        for z in range(nz):
            yfit = profile[slc, y, z]
            popt, pcov = curve_fit(fit_func, xfit, yfit, p0 = (nx//2, phi0, zeta))
            out[y, z] = popt[0]

    return out

def interface_height(profile, zero = False, level = 0, method = "direct", thermodynamic_params = None):
    """
    Calculates the interface height in the x-direction for a system, with an optional centering around zero.

    This function wraps the `ih_direct` function, providing the additional capability to center the interface height around zero. The interface height is determined for a flat interface across the yz plane.

    Parameters
    ----------
    profile : numpy.ndarray
        A 3D NumPy array representing the order parameter of the system. The array has dimensions `(nx, ny, nz)`, where `nx` is the size in the x-direction and `ny` and `nz` are the sizes in the y- and z-directions, respectively.
    zero : bool, optional
        If `True`, the output interface heights are centered around zero. Defaults to `False`.
    level: float, optional
        Sets the contour level for the find_contours function to calculate the interface height. Defaults to 0
    method: str, optional
        Sets which interface height calculation technique to use. "profile" uses a tanh fit while direct uses the find_contours method from skimage. Defaults to "direct"
    thermodynamic_params: list, optional
        Provides thermodynamic parameters for the "profile" fit to work. Is not used for "direct". Default is None.

    Returns
    -------
    numpy.ndarray
        A 2D NumPy array of shape `(ny, nz)` containing the interface heights in the yz plane. If `zero=True`, the heights are centered around zero.

    Notes
    -----
    - The centering is achieved by subtracting `(nx - 1) / 2` from the interface heights when `zero=True`.
    - The function uses `ih_direct` to calculate the raw interface heights.

    Examples
    --------
    Calculate the interface height without centering:

    >>> import numpy as np
    >>> profile = np.random.random((50, 30, 30)) - 0.5  # Example 3D order parameter
    >>> heights = interface_height(profile)
    >>> print(heights.shape)
    (30, 30)

    Calculate the interface height with centering around zero:

    >>> heights_centered = interface_height(profile, zero=True)
    >>> print(heights_centered.shape)
    (30, 30)
    """
    nx, ny, nz = profile.shape
    slc = np.s_[nx//4:3*nx//4, :, :]
    height_func = np.zeros((ny, nz))

    yraw = profile[slc].copy()
    zero_factor = (nx - 1)/2

    if method == "direct":
        height_func = ih_direct(yraw, level) + nx//4
    elif method == "profile":
        if thermodynamic_params is None:
            raise ValueError("Thermodynamic parameters not specified while using profile fit method")
        height_func = ih_profile_fit(profile, *thermodynamic_params)
    else:
        raise ValueError(f"method = {method} is invalid. Choose either direct or profile method options")
    
    height_func = height_func - zero_factor if zero else height_func
    
    return height_func

def calculate_surface_tension_profile(x, profile, kappa):
    grad_phi_x = np.gradient(profile)
    grad_square = np.power(grad_phi_x, 2)
    sigma_real = kappa*np.trapz(grad_square, x = x)
    return sigma_real