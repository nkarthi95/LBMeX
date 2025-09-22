import numpy as np
from scipy.optimize import curve_fit

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

def sharpen_profile(rho, rhoh=1, rhol=0, d=5/2):
    rhot = 0.5*(rhoh+rhol)
    return np.maximum(np.minimum((d*rho+(1-d)*rhot),rhoh),rhol)

def stokes_mobility(R, boxDims, alpha = 1, tau_r = 0.7886751345948129, rho = 1):
    L = boxDims.min()
    P = 1 - 2.84*R/L
    eta = rho*(1/3)*(tau_r - 0.5)
    fn = (6+4*alpha)/(1+alpha)
    return P/(fn*np.pi*eta*R)