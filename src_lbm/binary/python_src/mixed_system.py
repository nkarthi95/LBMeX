import numpy as np

def spherically_averaged_structure_factor(data, thermo_model, scale_factor = 1, func = None, shift = True, cs = True):

    """
    Calculates the spherically averaged structure factor for volumetric data.

    This function computes the spherically averaged structure factor from volumetric structure factor data. The function performs Fourier transformation and averaging over spherical shells in k-space, where the k-values are calculated based on the provided data. Additionally, the thermodynamic equation of state, scaling factor, and optional rescaling functions can be applied.

    Parameters
    ----------
    data : numpy.ndarray
        A 3D or 2D array containing the volumetric structure factor data. For 3D data, the function expects a 3D array, while for 2D data, a 2D array is accepted.

    thermo_model : object
        A custom thermodynamic model class instance that defines the thermodynamic equation of state, providing methods like `cs2k` or `mu_ck` for the k-dependent parameters. This model is used to adjust the structure factor based on thermodynamic properties.

    scale_factor : float, optional
        A factor by which the data is divided to scale the structure factor. The default value is `1`.

    func : callable, optional
        A lambda function for rescaling the data, which should take the data and the k-dependent values as input. The default value is `None`.

    shift : bool, optional
        A boolean that controls whether the 0-point of the Fourier transform is at the center (`True`) or the left (`False`). The default value is `True`.

    cs : bool, optional
        A boolean that determines whether the k-dependent parameter is calculated using `cs2k` (`True`) or `mu_ck` (`False`). The default value is `True`.

    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing two 1D arrays:
        - The first array contains the spherically averaged k-values.
        - The second array contains the spherically averaged structure factor values.

    Examples
    --------
    >>> data = np.random.random((64, 64, 64))  # Example 3D structure factor data
    >>> thermo_model = CustomThermoModel()  # Your custom thermodynamic model
    >>> k, S = spherically_averaged_structure_factor(data, thermo_model, scale_factor=1, shift=True, cs=True)

    Notes
    -----
    - The function calculates the spherically averaged structure factor by computing the k-values and structure factor values over spherical shells in k-space.
    - The `scale_factor` and `func` parameters allow for flexible scaling and rescaling of the structure factor.
    - The `shift` parameter controls whether the Fourier transform centers the zero-frequency component in the middle or at the left of the data.
    - The thermodynamic model (`thermo_model`) must provide functions like `cs2k` or `mu_ck`, which return k-dependent values needed for rescaling the structure factor.

    References
    ----------
    """
    L = min(data.shape)
    S = data.copy()

    if shift:
        freqs = np.fft.fftshift(np.fft.fftfreq(L))
    else:
        freqs = np.fft.fftfreq(L)
    if len(data.shape) == 3:
        kx, ky, kz = np.meshgrid(*tuple([2*np.pi*freqs for L in [L, L, L]]), indexing='ij')
        k = np.stack([kx, ky, kz], axis = -1)
    elif len(data.shape) == 2:
        kx, ky = np.meshgrid(*tuple([2*np.pi*freqs for L in [L, L]]), indexing='ij')
        k = np.stack([kx, ky], axis = -1)
     
    k1 = np.linalg.norm(k, axis=-1).flatten()
    
    if func is not None:
        if cs:
            S = func(S, thermo_model.cs2k(kx, ky, kz))
        else:
            S = func(S, thermo_model.mu_ck(kx, ky, kz))
    S /= scale_factor

    # test[slc] /= test.sum()
    # S[L//2, L//2, L//2] /= S.sum()

    S1 = S.flatten()
    kmin = 2*np.pi/L # sampling frequency
    where = np.s_[:]#np.where(k1<=kmax)
    bins = np.arange(L//2+1)*kmin # kmax+1 for bin_edges: len(bins)=len(hist)+1
    
    shells = np.histogram(k1[where], bins, weights=S1[where])[0]
    counts = np.histogram(k1[where], bins)[0]
    return (bins[:-1]+bins[1:])/2, shells/counts

def cart2sph(x,y,z):
    """
    Converts Cartesian coordinates to spherical coordinates.

    This function converts a 3D grid of Cartesian coordinates (x, y, z) into their corresponding spherical coordinates (radius, azimuth, and elevation). The conversion follows standard spherical coordinate transformations:
    - The radius is the distance from the origin.
    - The azimuth is the angle in the xy-plane, relative to the x-axis.
    - The elevation is the angle from the xy-plane, relative to the z-axis.

    Parameters
    ----------
    x : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the x-coordinates of the Cartesian grid.

    y : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the y-coordinates of the Cartesian grid.

    z : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the z-coordinates of the Cartesian grid.

    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing three 3D numpy arrays:
        - The first array is the radial distance (r) from the origin at each grid point.
        - The second array is the azimuthal angle (azimuth), the angle in the xy-plane from the x-axis.
        - The third array is the elevation angle (elevation), the angle from the xy-plane to the z-axis.

    Examples
    --------
    >>> x = np.random.random((10, 10, 10))  # Example 3D x-coordinate grid
    >>> y = np.random.random((10, 10, 10))  # Example 3D y-coordinate grid
    >>> z = np.random.random((10, 10, 10))  # Example 3D z-coordinate grid
    >>> r, azimuth, elevation = cart2sph(x, y, z)

    Notes
    -----
    - The `azimuth` is calculated using `arctan2(y, x)`, which gives the angle in the xy-plane, taking into account the quadrant of (x, y).
    - The `elevation` is calculated using `arctan2(z, sqrt(x^2 + y^2))`, giving the angle relative to the xy-plane.
    - The radial distance `r` is computed as `sqrt(x^2 + y^2 + z^2)`.

    """
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    return r, azimuth, elevation

def sph2cart(azimuth,elevation,r):
    """
    Converts spherical coordinates to Cartesian coordinates.

    This function converts spherical coordinates (radius, azimuth, elevation) into their corresponding Cartesian coordinates (x, y, z). The conversion follows standard spherical coordinate transformations:
    - The radius is the distance from the origin.
    - The azimuth is the angle in the xy-plane, relative to the x-axis.
    - The elevation is the angle from the xy-plane, relative to the z-axis.

    Parameters
    ----------
    azimuth : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the azimuthal angles (in radians) of the spherical grid.

    elevation : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the elevation angles (in radians) of the spherical grid.

    r : numpy.ndarray
        A 3D numpy array of shape (nx, ny, nz) containing the radial distances (r) from the origin in the spherical grid.

    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing three 3D numpy arrays:
        - The first array is the x-coordinate in the Cartesian grid.
        - The second array is the y-coordinate in the Cartesian grid.
        - The third array is the z-coordinate in the Cartesian grid.

    Examples
    --------
    >>> azimuth = np.random.random((10, 10, 10))  # Example 3D azimuthal angle grid (in radians)
    >>> elevation = np.random.random((10, 10, 10))  # Example 3D elevation angle grid (in radians)
    >>> r = np.random.random((10, 10, 10))  # Example 3D radial distance grid
    >>> x, y, z = sph2cart(azimuth, elevation, r)

    Notes
    -----
    - The azimuth is assumed to be the angle in the xy-plane, and the elevation is the angle from the xy-plane to the z-axis.
    - The conversion uses the following equations:
    - `x = r * cos(elevation) * cos(azimuth)`
    - `y = r * cos(elevation) * sin(azimuth)`
    - `z = r * sin(elevation)`
    """
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

def radial_equilibration(data, thermo_model, radius = 1, scale_factor = 1, func = None, cs = True):
    """
    Calculates the radial average of the input data at the same radius.

    This function computes the average of the input data at a specified radius in spherical coordinates. The data is first transformed from Cartesian to spherical coordinates. If a post-processing function is provided, it is applied to the data. The function allows for scaling using either the speed of sound or chemical potential, defined by the `cs` parameter.

    Parameters
    ----------
    data : numpy.ndarray
        A 3D numpy array representing the input data on a grid. The data is assumed to be in Cartesian coordinates.

    thermo_model : object
        A custom class instance that specifies the thermodynamic equation of state. This model is used to calculate scaling factors based on either the speed of sound or the chemical potential.

    radius : float, optional
        A scalar value defining the radius over which to compute the average. Defaults to 1.0.

    scale_factor : float, optional
        A scaling factor to divide all the data. Defaults to 1.

    func : lambda function, optional
        A lambda function that performs post-processing on the data. It is commonly used to apply a k-dependent scaling. Defaults to None.

    cs : bool, optional
        A boolean that determines which scaling factor to use:
        - If True, the speed of sound is used as the scaling factor.
        - If False, the chemical potential is used as the scaling factor. Defaults to True.

    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing three 1D numpy arrays:
        - The first array contains the azimuthal angles (t).
        - The second array contains the polar angles (p).
        - The third array contains the averaged data values at the corresponding points.

    Examples
    --------
    >>> data = np.random.random((100, 100, 100))  # Example 3D data array
    >>> thermo_model = MyThermodynamicModel()  # Example thermodynamic model class
    >>> t, p, out = radial_equilibration(data, thermo_model, radius=1, scale_factor=1, func=None, cs=True)

    Notes
    -----
    - The input data is assumed to be in Cartesian coordinates and is converted to spherical coordinates during the process.
    - The azimuthal (t) and polar (p) angles are returned along with the averaged data values at the specified radius.
    - The `func` argument allows for post-processing of the data, such as applying a k-dependent scaling factor, which is common in Fourier space calculations.
    - The `radius` parameter is used to select the specific spherical shell for averaging, and the default value of 1 corresponds to the unit sphere.
    """

    S = data.copy()
    L = min(S.shape)
    freqs = np.fft.fftshift(np.fft.fftfreq(L))
    kx, ky, kz = np.meshgrid(*tuple([2*np.pi*freqs for L in [L, L, L]]), indexing='ij')

    r, t, p = cart2sph(kx, ky, kz)

    if func is not None:
        if cs:
            S = func(S, thermo_model.cs2k(kx, ky, kz))
        else:
            S = func(S, thermo_model.mu_ck(kx, ky, kz))
    S /= scale_factor
    # S[L//2, L//2, L//2] /= S.sum()
    
    idxs = np.isclose(r, radius, atol = 2*np.pi/L)
    t = t[idxs]
    p = p[idxs]
    out = S[idxs]

    return t, p, out

def make_bins(to_bin1, binsize, to_bin2 = None):
    """
    Performs binning of input arrays to a specified bin size.

    This function takes two 1D numpy arrays and bins them into specified bins. The first input array, `to_bin1`, is binned into a specified number of bins, and if a second array, `to_bin2`, is provided, it will be binned to the same size as the first array. The function returns one or two binned arrays depending on whether `to_bin2` is provided.

    Parameters
    ----------
    to_bin1 : numpy.ndarray
        A 1D numpy array of data to be binned. The array will be binned based on the specified `binsize`.

    binsize : int
        The number of bins to divide the data into. The resulting binned arrays will have this size.

    to_bin2 : numpy.ndarray, optional
        A second 1D numpy array of data to be binned to the same size as `binsize`. This parameter is optional. If provided, it will be binned and returned alongside the first binned array. Defaults to None.

    Returns
    -------
    tuple
        If `to_bin2` is not specified, returns a tuple containing:
        - A 1D numpy array of the bin edges.
        - A 1D numpy array of the binned values of `to_bin1`.

        If `to_bin2` is specified, returns:
        - A 1D numpy array of the binned values of `to_bin1`.
        - A 1D numpy array of the binned values of `to_bin2`.

    Examples
    --------
    >>> to_bin1 = np.random.random(1000)  # Example data
    >>> to_bin2 = np.random.random(1000)  # Example second data
    >>> bins, binned_data = make_bins(to_bin1, 50)  # Binning to_bin1 into 50 bins
    >>> binned_data1, binned_data2 = make_bins(to_bin1, 50, to_bin2)  # Binning both arrays

    Notes
    -----
    - The bins are created using `numpy.linspace` to define the bin edges between the minimum and maximum of `to_bin1`.
    - The function uses `numpy.digitize` to assign each data point in `to_bin1` (and `to_bin2`, if provided) to a bin. The `np.add.at` function is then used to sum the data points in each bin, followed by averaging the values in each bin.
    - If `to_bin2` is provided, both `to_bin1` and `to_bin2` are binned and returned as separate arrays.
    """


    bins = np.linspace(to_bin1.min(), to_bin1.max(), binsize)

    out1 = np.zeros(binsize)
    shell = np.digitize(to_bin1, bins = bins, right = True)
    np.add.at(out1, shell, to_bin1)
    unique, counts = np.unique(shell, return_counts=True)
    out1 = out1[unique]
    out1 /= counts

    if to_bin2 is None:
        return bins, out1
    else:
        out2 = np.zeros(binsize)
        np.add.at(out2, shell, to_bin2)
        unique, counts = np.unique(shell, return_counts=True)
        out2 = out2[unique]
        out2 /= counts
        return out1, out2