import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import center_of_mass

def pressure_jump(pressure):
    """
    Calculates the pressure jump between the interior of a droplet and the corners of the simulation domain.

    The function computes the pressure difference (`dP`) by averaging the pressure in the interior (a small slice at the center of the domain) and the pressure at the corners of the domain. This is useful for analyzing the pressure distribution in systems such as droplets or bubbles.

    Parameters
    ----------
    pressure : numpy.ndarray
        A 3D NumPy array representing the scalar pressure field of the system. The array has dimensions `(nx, ny, nz)`.

    Returns
    -------
    float
        The pressure jump (`dP`) between the interior of the droplet and the corners of the simulation domain.

    Notes
    -----
    - The interior pressure is calculated as the average of a 3x3x3 slice centered at the middle of the domain.
    - The exterior pressure is calculated as the average of all values at the corners of the simulation domain.

    Examples
    --------
    Compute the pressure jump for a given pressure field:

    >>> import numpy as np
    >>> pressure = np.random.random((50, 50, 50))  # Example 3D pressure field
    >>> dP = pressure_jump(pressure)
    >>> print(dP)
    0.023456789
    """

    nx, ny, nz = pressure.shape
    center_slc = np.s_[nx//2-1:nx//2+2, ny//2-1:ny//2+2, nz//2-1:nz//2+2]
    edge_slc = np.s_[0:nx:nx-1, 0:ny:ny-1, 0:nz:nz-1]
    dP = pressure[center_slc].mean() - pressure[edge_slc].mean()
    return dP

def droplet_mass(OutArray):
    """
    Calculates the mass of a droplet.

    This function calculates the mass of a droplet by summing the values of the order parameter encoded in a 3D NumPy array. The result corresponds to the total mass of the droplet, assuming the mass is proportional to the order parameter.

    Parameters
    ----------
    OutArray : numpy.ndarray
        A 3D NumPy array representing the order parameter of the system, where the values correspond to the concentration of the droplet at each point.

    Returns
    -------
    float
        The mass of the droplet, calculated as the sum of the values in the `OutArray`.

    Examples
    --------
    Calculate the mass of a droplet given a 3D order parameter array:

    >>> import numpy as np
    >>> OutArray = np.random.random((50, 50, 50))  # Example 3D order parameter array
    >>> mass = droplet_mass(OutArray)
    >>> print(mass)
    1234.56  # Example mass value
    """
    sum = np.sum(OutArray)
    return sum

def droplet_radius_mass(density, Vp = 0, np_sphere = 0, rho_sphere = 1):
    """
    Calculate the radius of a droplet from the density field and particle parameters.

    This function calculates the radius of a droplet based on the density distribution of two 
    density fields. It also considers particle volume, number of particles, and the density of 
    the particles in the calculation. The radius is computed using the mass and density difference 
    within the system.

    :param density: 
        A 3D numpy array representing the density field of the system, where each value encodes 
        the local density at a given point in the simulation box.
    :type density: numpy.ndarray

    :param Vp: 
        The volume of a single particle. This parameter is used to estimate the mass of the droplet 
        based on the number of particles in the system. Default is 0.
    :type Vp: float, optional

    :param np_sphere: 
        The total number of particles in the droplet or the system. This is used to calculate the 
        mass contribution of the particles in the droplet. Default is 0.
    :type np_sphere: int, optional

    :param rho_sphere: 
        The density of the particles. This is used to calculate the total mass of the particles 
        contributing to the droplet. Default is 1.
    :type rho_sphere: float, optional

    :return: 
        The radius of the droplet, calculated based on the mass and the density difference 
        between the droplet and the surrounding medium.
    :rtype: float

    :note: 
        - The function assumes that the density field is centered on the droplet, and it calculates 
          the radius by considering the mass of the droplet and its density contrast with the surrounding 
          medium.
        - If `density` is an integer, the function returns `NaN`, as this indicates invalid input.

    :example:
        >>> density = np.random.random((10, 10, 10))  # Example density field
        >>> droplet_radius(density, Vp=1, np_sphere=100, rho_sphere=1.5)
        1.25  # Example output for droplet radius
    """
    if isinstance(density, int):
        return np.nan    
    else:
        nx, ny, nz = density.shape
        center_slc = np.s_[nx//2-1:nx//2+2, ny//2-1:ny//2+2, nz//2-1:nz//2+2]
        edge_slc = np.s_[0:nx:nx-1, 0:ny:ny-1, 0:nz:nz-1]
        # center = tuple([ l//2 for l in density.shape ])
        
        rho_d = density[center_slc].mean()
        rho_m = density[edge_slc].mean()
        # mass = np.sum(density - rho_m) + 0.5*Vp*np_sphere*rho_sphere
        mass = droplet_mass(density - rho_m) + 0.5*Vp*np_sphere*rho_sphere
        R = (3./4./np.pi*mass/(rho_d-rho_m))**(1./3.)
        return R
    
# def droplet_radius_profile_1d(density, center=None, fit='tanh', bins=100, binned = False):
#     """
#     Estimate droplet radius and interface width from a 3D density field
#     by fitting a 1D radial density profile to a tanh or erf function.
    
#     Parameters
#     ----------
#     density : ndarray
#         3D array of concentration or density values.
#     center : tuple or None
#         Optional center of mass (x, y, z). If None, center is computed from density.
#     fit : str
#         'tanh' or 'erf' for the fitting function.
#     bins : int
#         Number of radial bins for the profile.

#     Returns
#     -------
#     popt : array
#         Fit parameters: [radius R, amplitude p0, baseline c, interfacial width xi]
#     r_bin_centers : array
#         Radii of the profile bins (for plotting, optional)
#     rho_r : array
#         Radial density profile used in the fit
#     """

#     X, Y, Z = np.indices(density.shape)

#     # Estimate center if not provided
#     if center is None:
#         cmx, cmy, cmz = center_of_mass(density)
#     else:
#         cmx, cmy, cmz = center

#     # Compute radial distances from center
#     r = np.sqrt((X - cmx)**2 + (Y - cmy)**2 + (Z - cmz)**2).flatten()
#     rho = density.flatten()

#     if binned:
#         # Bin by radius
#         r_bins = np.linspace(0, r.max(), bins + 1)
#         r_bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
#         rho_r = np.zeros(bins)

#         for i in range(bins):
#             mask = (r >= r_bins[i]) & (r < r_bins[i+1])
#             rho_r[i] = np.mean(rho[mask]) if np.any(mask) else np.nan

#         # Remove NaNs (can happen if outer bins are empty)
#         valid = ~np.isnan(rho_r)
#         r_bin_centers = r_bin_centers[valid]
#         rho_r = rho_r[valid]
#     else:
#         idxs_ascending = np.argsort(r)
#         r_bin_centers = r[idxs_ascending]
#         rho_r = rho[idxs_ascending]

#     # Define fitting function
#     if fit == 'tanh':
#         def fit_func(r, R, p0, c, xi):
#             return p0 * np.tanh((R - r) / (np.sqrt(2) * xi)) + c
#     elif fit == 'erf':
#         from scipy.special import erf
#         def fit_func(r, R, p0, c, xi):
#             return p0 * erf((R - r) / (np.sqrt(2) * xi)) + c
#     else:
#         raise ValueError("Invalid fit type. Use 'tanh' or 'erf'.")

#     # Guess initial parameters
#     p0_guess = 0.5 * (np.max(rho_r) - np.min(rho_r))
#     c_guess = 0.5 * (np.max(rho_r) + np.min(rho_r))
#     R_guess = r_bin_centers[np.argmin(np.abs(rho_r - c_guess))]
#     xi_guess = 0.3
#     guess = (R_guess, p0_guess, c_guess, xi_guess)

#     # Fit
#     popt, _ = curve_fit(fit_func, r_bin_centers, rho_r, p0=guess)

#     return popt, r_bin_centers, rho_r

def get_box(rho):
    #box = np.stack(np.meshgrid(*[range(L) for L in rho.shape], indexing='ij'), axis=-1)-np.asarray(rho.shape)/2+0.5
    box = np.moveaxis(np.indices(rho.shape), 0, -1) - np.asarray(rho.shape)/2 + 0.5
    return box

def droplet_profile(rho):
    profile = lambda r, R, alpha, rhoh, rhol: rhol+0.5*(rhoh-rhol)*(1+np.tanh((R-np.sqrt(r**2))/(alpha/2)))
    rs = np.linalg.norm(get_box(rho), axis=-1)
    p = curve_fit(profile, rs.flatten(), rho.flatten())[0]
    #R, alpha, rhoh, rhol = p
    return p

def gyration_tensor(OutArray, cm):
    """
    Calculate the gyration tensor of a 3D array with respect to its center of mass.

    This function computes the gyration tensor of a 3D numpy array `OutArray` using the provided 
    center of mass `cm`. The gyration tensor provides a measure of the spatial distribution of mass 
    around the center of mass.

    :param cm: 
        A 1D array or list of size 3 representing the center of mass of the array.
    :type cm: array-like

    :param OutArray: 
        A 3D numpy array where each element represents a scalar mass density at that position.
    :type OutArray: numpy.ndarray

    :return: 
        A (3, 3) numpy matrix representing the gyration tensor of `OutArray`.
    :rtype: numpy.ndarray

    :note: 
        - The gyration tensor is normalized by the total mass (sum of `OutArray`) and is calculated 
          using the second moment of the position vectors relative to the center of mass.
        - This tensor is useful for characterizing the shape and size of spatial distributions.

    :example:
        >>> cm = np.array([5.0, 5.0, 5.0])  # Center of mass
        >>> OutArray = np.random.random((10, 10, 10))  # Example density field
        >>> S = gyration_tensor(cm, OutArray)
        >>> print(S)
        [[0.33, 0.01, 0.02],
         [0.01, 0.35, 0.03],
         [0.02, 0.03, 0.37]]  # Example output gyration tensor
    """
    ind = np.transpose(np.indices(OutArray.shape), axes=(1,2,3,0))
    pos = ind - cm
    rr = np.einsum('...m,...n->...mn',pos,pos)
    S = np.einsum('ijk,ijk...',OutArray,rr)/np.sum(OutArray)
    return S

def gyration_tensor(rho, cm=(0,0,0)):
    box = np.array(rho.shape)
    r = get_box(rho) - cm
    r -= box*(r/box+np.sign(r)*0.5).astype(int) # minimum image convention
    S = np.einsum('ijk,ijka,ijkb', rho, r, r)/np.sum(rho)
    return S

def center_of_mass_us(rho):
    com = np.einsum('ijk,ijka', rho, get_box(rho))/np.sum(rho)
    return com

def axial_radii(field, cm = None):
    """
    Calculates the fluctuations in the principal radii of a droplet.

    This function computes the fluctuations in the principal radii of a droplet by calculating the gyration tensor, extracting its eigenvalues, and then computing the variation in each of the droplet's principal axes (x, y, and z). The fluctuations are defined relative to the ideal case of no fluctuation (where the fluctuation value is 1). The output is a 3D vector representing the fluctuations in the x, y, and z axes.

    Parameters
    ----------
    field : numpy.ndarray
        A 3D NumPy array representing the order parameter of the system, encoding the droplet shape and composition. The array is used to compute the center of mass and gyration tensor.

    Returns
    -------
    numpy.ndarray
        A 1D NumPy array of size 3, representing the fluctuation in the droplet's radii along the x, y, and z directions, respectively. A value of 1 corresponds to no fluctuation, and values greater than 1 indicate increased fluctuation.

    Examples
    --------
    Calculate the fluctuations in the principal radii of a droplet:

    >>> import numpy as np
    >>> field = np.random.random((50, 50, 50))  # Example 3D order parameter field
    >>> fluctuations = axial_radii(field)
    >>> print(fluctuations)
    [1.05 0.98 1.12]  # Example fluctuations in x, y, and z directions

    Notes
    -----
    - The fluctuations are calculated by first obtaining the gyration tensor and its eigenvalues.
    - The fluctuation in each direction is then calculated based on the ratio of the first eigenvalue to the geometric mean of the other two eigenvalues.
    - This function assumes the input field represents the entire droplet system, and that the gyration tensor can be computed from it.

    """
    if cm is None:
        cm = center_of_mass_us(field)
    
    gr = gyration_tensor(field, cm) # Calculating the gyration tensor of the droplet
    egr = np.linalg.eigvals(gr) # calculating the unordered eigenvalues of the gyration tensor
    # egr = np.sqrt(egr) # calculating the unordered eigenvalues of the gyration tensor
    da = np.power(egr[0], 1/3)/np.power(np.prod(egr[[1,2]]), 1/6) # calculating the variations in dx
    db = np.power(egr[1], 1/3)/np.power(np.prod(egr[[0,2]]), 1/6) # calculating the variations in dy
    dc = np.power(egr[2], 1/3)/np.power(np.prod(egr[[0,1]]), 1/6) # calculating the variations in dz  
    return np.array([da, db, dc])

def droplet_fluctuations(fluctuations, temp = 1e-7):
    """
    Calculates the surface tension of a droplet from fluctuations in its principal radii.

    This function calculates the surface tension of a droplet based on the fluctuations in its principal radii along the x, y, and z directions. The fluctuations are used to compute surface tension from two distinct harmonic modes, which are characterized by sums and differences of the radii fluctuations. The function returns the surface tension for each mode as a list of two values.

    Parameters
    ----------
    fluctuations : numpy.ndarray
        A 2D NumPy array of shape (N, 3), where N is the number of observations, and 3 represents the principal radii fluctuations in the x, y, and z directions.
    temp : float, optional
        A thermalization parameter representing the temperature of the system. The default value is `1e-7`.

    Returns
    -------
    list of float
        A list of length 2, where the first element corresponds to the surface tension from the y20 harmonic mode, and the second element corresponds to the surface tension from the y22 harmonic mode.

    Examples
    --------
    Calculate the surface tension from droplet radius fluctuations:

    >>> fluctuations = np.random.random((100, 3))  # Example fluctuations (100 observations)
    >>> temp = 1.0  # Example temperature
    >>> surface_tension = droplet_fluctuations(fluctuations, temp)
    >>> print(surface_tension)
    [0.025 0.042]  # Example surface tension values for the y20 and y22 modes

    Notes
    -----
    - The surface tension for the y20 and y22 harmonic modes is calculated using the fluctuations of the principal radii along the x, y, and z axes.
    - The temperature parameter `temp` is used to scale the surface tension values and defaults to a very small value (`1e-7`) if not provided.
    - The final surface tension values are computed by averaging the squared sums and differences of the fluctuations in the respective directions.

    """
    sums = 0
    difs = 0

    for i in range(0, 2):
        for j in range(i+1, 3):
            sums += np.mean(np.power(fluctuations[:, i] + fluctuations[:, j], 2))
            difs += np.mean(np.power(fluctuations[:, i] - fluctuations[:, j], 2))

    sums *= 1/3
    difs *= 1/3

    y20 = 5*temp/(16*np.pi*sums)
    y22 = 15*temp/(16*np.pi*difs)

    return [y20, y22]

# def sharpen_droplet_interface(profile, dcf = 1/5, int_height = 0, max_val = 1, min_val = 0):
#     """
#     Applies an interface sharpening transformation to an n-dimensional numpy array.

#     This function modifies the input array `profile` by adjusting the interface 
#     sharpness based on the `dcf` parameter. The sharpening process ensures that 
#     values are constrained within a specified range, effectively enhancing the 
#     contrast at the interface.

#     Parameters
#     ----------
#     profile : numpy.ndarray
#         An n-dimensional numpy array representing the scalar field to be sharpened.
#     dcf : float, optional
#         A parameter that controls the strength of the interface sharpening. 
#         Default is 1/5.
#     int_height: float, optional
#         A parameter that identifies the crossover point defining the interface in a phase separated
#         system. Defaults to 0.

#     Returns
#     -------
#     numpy.ndarray
#         An n-dimensional numpy array after the interface sharpening transformation.

#     Notes
#     -----
#     The transformation is computed as:
#         density_filt = (profile + dcf - 0.5) / (2 * dcf)
#     Follows the method detailed in Equation 9 of https://doi.org/10.1063/5.0249847
#     """
#     density_filt = (profile + dcf - int_height)/(2*dcf)
#     density_filt = np.where(density_filt < max_val, density_filt, max_val)
#     density_filt = np.where(density_filt > min_val, density_filt, min_val)
#     return density_filt

def sharpen_profile(rho, rhoh=1, rhol=0, d=5/2):
    rhot = 0.5*(rhoh+rhol)
    return np.maximum(np.minimum((d*rho+(1-d)*rhot),rhoh),rhol)