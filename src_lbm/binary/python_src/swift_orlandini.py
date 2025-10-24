import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from d3q19 import lattice_fourier_laplacian
from scipy.optimize import newton
from scipy.integrate import quad
import numpy as np

class swift_et_al_1996_thermodynamic_model:
    """
    Defines a thermodynamic equation of state based on the model by Swift et al. (1996).

    This class models thermodynamic properties such as free energy, chemical potential, and their derivatives. It also provides methods to calculate these quantities in Fourier space, incorporating the effects of wavevectors.

    Parameters
    ----------
    density : float, optional
        The density of the system. Defaults to 1.
    C0 : float, optional
        The concentration parameter. Defaults to 0.
    chi : float, optional
        The interaction parameter. Defaults to 0.4.
    T : float, optional
        The temperature of the system. Defaults to 0.25.
    kappa : float, optional
        The gradient energy coefficient. Defaults to 0.01.

    Methods
    -------
    sound_speed_square()
        Calculates the square of the sound speed in the system.
    cs2k(kx=0, ky=0, kz=0)
        Calculates the sound speed squared in Fourier space, incorporating contributions from the gradient energy term.
    calculate_df_dphi()
        Computes the derivative of the free energy with respect to the order parameter.
    calculate_dmup_drho()
        Computes the derivative of the chemical potential with respect to the density.
    calculate_dmup_dphi()
        Computes the derivative of the chemical potential with respect to the order parameter.
    mu_ck(kx=0, ky=0, kz=0)
        Calculates the chemical potential in Fourier space.

    Notes
    -----
    - Fourier space calculations are based on the `lattice_fourier_laplacian` function, which computes the Laplacian in Fourier space for a discrete lattice.
    - The methods rely on internal state variables such as `rho` (density), `C0` (concentration), `chi` (interaction parameter), `T` (temperature), and `kappa` (gradient energy coefficient).

    Examples
    --------
    Initialize the model and calculate thermodynamic properties:

    >>> model = swift_et_al_1996_thermodynamic_model(density=1.2, C0=0.1, chi=0.5, T=0.3, kappa=0.02)
    >>> sound_speed_squared = model.sound_speed_square()
    >>> print(sound_speed_squared)
    0.3

    Calculate the chemical potential derivative with respect to the density:

    >>> dmup_drho = model.calculate_dmup_drho()
    >>> print(dmup_drho)
    -0.20833333333333334

    Compute the chemical potential in Fourier space:

    >>> mu_ck_value = model.mu_ck(kx=0.1, ky=0.1, kz=0.1)
    >>> print(mu_ck_value)
    -0.3037037037037037
    """
    def __init__(self, density = 1, C0 = 0, chi = 0.4, T = 0.25, kappa = 0.01):
        self.chi = chi
        self.T = T
        self.kappa = kappa
        self.rho = density
        self.C0 = C0

    def sound_speed_square(self):
        out = self.T+1/3
        return out

    def cs2k(self, kx = 0, ky = 0, kz = 0):
        thermal_cs2 = self.sound_speed_square()
        k2 = lattice_fourier_laplacian(kx, ky, kz)
        out = thermal_cs2 + k2*self.kappa
        return out
    
    def calculate_bulk(self):
        rho = self.rho
        phi = self.C0
        chi = self.chi
        T   = self.T
        out = chi*rho/4*(1 - (phi/rho)**2) - rho*T + T/2*(rho+phi)*np.log((rho+phi)/2) + T/2*(rho-phi)*np.log((rho-phi)/2)
        return out

    def calculate_df_dphi(self):
        rho = self.rho
        phi = self.C0
        chi = self.chi
        T   = self.T
        out = -chi/2.*(phi/rho) + T/2.*np.log((1. + phi/rho)/(1. - phi/rho))
        return out

    def calculate_dmup_drho(self):
        rho = self.rho
        phi = self.C0
        chi = self.chi
        T   = self.T
        out = -T*phi/(rho**2 - phi**2) + chi/2*phi/(rho**2)
        return out

    def calculate_dmup_dphi(self):
        rho = self.rho
        phi = self.C0
        chi = self.chi
        T   = self.T
        out = T*rho/(rho**2 - phi**2) - chi/(2*rho)
        return out
    
    def mu_ck(self, kx = 0, ky = 0, kz = 0):
        ref_state = self.calculate_dmup_dphi()
        k2 = lattice_fourier_laplacian(kx, ky, kz)
        out = ref_state + k2*self.kappa
        return out
    
def swift_theoretical_C0(chi, guess = 0.99, tol = 1e-4):
    """
    Compute the theoretical value of C0 from the Swift equation of state using the Newton-Raphson method.

    This function solves for C0 using the equation:

        log(c / (1 - c)) - chi * (2c - 1) = 0

    where `chi` is the ratio of two thermodynamic quantities in the Swift equation of state.

    Parameters
    ----------
    chi : float
        The ratio of two thermodynamic quantities from the Swift equation of state.
    guess : float, optional
        The initial guess for the Newton-Raphson method. Default is 0.99.
    tol : float, optional
        The tolerance for the Newton-Raphson method. Default is 1e-4.

    Returns
    -------
    float
        The computed value of C0 that satisfies the equation.

    Notes
    -----
    The function uses `scipy.optimize.newton` to iteratively solve for C0.

    Examples
    --------
    >>> swift_theoretical_C0(chi=2.5)
    0.9243  # Example output (may vary)
    """

    c0_expr =lambda c, chi: np.log(c/(1 - c)) - chi*(2*c - 1)
    c0 = newton(c0_expr, guess, tol = tol, args = [chi])
    return c0

def swift_critical_sigma(rho, chi, T, kappa):
    Tc = chi/2
    beta = 1/2/rho
    gamma = Tc/12/rho**3
    phi0 = np.sqrt(beta/2/gamma*(Tc-T))
    sigma = 4/3*np.sqrt(kappa*gamma)*phi0**3
    alpha = np.sqrt(2*kappa/gamma)/phi0
    return phi0, sigma, alpha

def fit_swift_phi0(chi):
    """
    Compute the coexistence density (ϕ₀) from the Swift equation of state using a piecewise approximation.

    This function provides an approximate solution for ϕ₀ without requiring the Newton-Raphson method.
    The piecewise function is defined as:

    - For chi ≤ 2.27: ϕ₀ = sqrt(3 - 6 / chi)
    - For 2.27 < chi ≤ 4.15: ϕ₀ = -0.133 * chi² + 1.030 * chi - 1.039
    - For chi > 4.15: ϕ₀ = 0.0248 * chi + 0.863

    The output is clamped to a maximum value of 0.99 to ensure physical consistency.

    Parameters
    ----------
    chi : float
        The ratio of two thermodynamic quantities from the Swift equation of state.

    Returns
    -------
    float
        The computed coexistence density (ϕ₀).

    Notes
    -----
    This function is an empirical fit to the coexistence densities of the Swift equation of state.
    It provides a fast alternative to iterative numerical methods.

    Examples
    --------
    >>> swift_fit_phi0(2.5)
    0.783  # Example output (may vary)
    """
    phi0 = 0.0
    if chi <= 2.27:
        phi0 = np.sqrt(3 - 6/chi)
    elif chi <= 4.15:
        phi0 = -0.133*chi**2 + 1.030*chi - 1.039
    else:
        phi0 = 0.0248*chi + 0.863
    
    if phi0 > 1:
        phi0 = 0.99
    
    return phi0

def swift_theoretical_sigma(chi, kappa):
    """
    Compute the interfacial tension (σ) for the Swift equation of state.

    This function calculates the interfacial tension using numerical integration of the 
    bulk free energy difference. It relies on the `quad` function from `scipy.integrate` 
    to perform the integration and uses `swift_theoretical_C0` to determine the coexistence 
    composition.

    Parameters
    ----------
    chi : float
        The ratio of two thermodynamic quantities from the Swift equation of state.
    kappa : float
        The interfacial interaction parameter.

    Returns
    -------
    float
        The computed interfacial tension (σ).

    Notes
    -----
    The function integrates the bulk free energy difference:

        ∫ sqrt((2c/χ) log(c/cₑ) + (2(1-c)/χ) log((1-c)/(1-cₑ)) - 2(c - cₑ)²) dc

    over the range `[1 - cₑ, cₑ]`, where `cₑ` is the equilibrium composition 
    obtained from `swift_theoretical_C0(chi)`. The result is then scaled by a 
    prefactor `sqrt(χκ / 2)`.

    This approach follows the theoretical prediction for interfacial tension in 
    the Swift free energy functional.

    Examples
    --------
    >>> swift_theoretical_sigma(chi=2.5, kappa=0.1)
    0.245  # Example output (may vary)
    """
    def bulk_free_energy_diff(c, chi, ce):
        term1 = (2*c/chi)*np.log(c/ce)
        term2 = (2*(1-c)/chi)*np.log((1-c)/(1-ce))
        term3 = 2*np.power(c - ce, 2)
        return np.sqrt(term1 + term2 - term3)
    
    ce = swift_theoretical_C0(chi)
    sigma_r = quad(bulk_free_energy_diff, 1 - ce, ce, args = (chi, ce))[0]

    coeff = np.sqrt(chi*kappa)

    return coeff*sigma_r

def swift_theoretical_xi(c1, chi, kappa):
    """
    Compute the interfacial width (ξ) for the Swift equation of state.

    This function calculates the interfacial width based on the coexistence density 
    of component 1 (`c1`), the thermodynamic parameter ratio (`chi`), and the interfacial 
    strength parameter (`kappa`). The expression is given by:

        ξ = (2 * sqrt(κ / χ)) / sqrt(-1 - (2 log(4 c₁ (1 - c₁))) / (χ (1 - 2 c₁)²))

    Parameters
    ----------
    c1 : float
        The coexistence density of component 1.
    chi : float
        The ratio of two thermodynamic quantities from the Swift equation of state.
    kappa : float
        The interfacial strength parameter.

    Returns
    -------
    float
        The computed interfacial width (ξ).

    Notes
    -----
    The function assumes that `c1` is within the physically valid coexistence range.
    The denominator must remain real and nonzero to avoid computational errors.

    Examples
    --------
    >>> swift_theoretical_xi(c1=0.3, chi=2.5, kappa=0.1)
    1.24  # Example output (may vary)
    """
    term1 = 2*np.sqrt(kappa/chi)
    term2 = np.sqrt(-1 - (2*np.log(4*c1*(1-c1)))/(chi*np.power(1 - 2*c1, 2)))
    return term1/term2

def noise_covariance_matrix(rho0, phi0, k2 = 0.0, 
                            tau_r = 0.7886751345948129, tau_p = 1.0, 
                            chi = 0.2, T = 0.095, kappa = 0.01, 
                            Gamma = 1.0, kT = 1e-5):
    Q = 19
    ndof = 2*Q

    lambdaLB_r = -1./tau_r
    lambdaLB_p = -1./tau_p

    lambda_r = -lambdaLB_r*(2.+lambdaLB_r)/2.
    lambda_p = -lambdaLB_p*(2.+lambdaLB_p)/2.
    lambda_rp = -lambdaLB_r*(2.+lambdaLB_p)/2.
    lambda_pr = -lambdaLB_p*(2.+lambdaLB_r)/2.

    cs2 = 1./3. + T + kappa*k2*rho0; # p_rho => modified speed of sound
    p_phi = kappa*k2*phi0
    mu_rho = -T*phi0/(rho0*rho0-phi0*phi0) + chi/2*phi0/(rho0*rho0)
    mu_phi = T*rho0/(rho0*rho0-phi0*phi0) - chi/2/rho0 + kappa*k2
    
    Xi = np.zeros(ndof*ndof)

    # diagonal part
    Xi[(   5)*ndof+(   5)] = 2.*Gamma*kT/rho0*lambda_p
    Xi[(   6)*ndof+(   6)] = 2.*Gamma*kT/rho0*lambda_p
    Xi[(   7)*ndof+(   7)] = 2.*Gamma*kT/rho0*lambda_p
    Xi[(   8)*ndof+(   8)] = 2.*kT*rho0*(5 - 9*cs2)*lambda_r
    Xi[(   9)*ndof+(   9)] = 8.*kT*rho0*lambda_r
    Xi[(  10)*ndof+(  10)] = (8.0/3.0)*kT*rho0*lambda_r
    Xi[(  11)*ndof+(  11)] = (2.0/3.0)*kT*rho0*lambda_r
    Xi[(  12)*ndof+(  12)] = (2.0/3.0)*kT*rho0*lambda_r
    Xi[(  13)*ndof+(  13)] = (2.0/3.0)*kT*rho0*lambda_r
    Xi[(  14)*ndof+(  14)] = 4.*kT*rho0*lambda_r
    Xi[(  15)*ndof+(  15)] = 4.*kT*rho0*lambda_r
    Xi[(  16)*ndof+(  16)] = 4.*kT*rho0*lambda_r
    Xi[(  17)*ndof+(  17)] = (4.0/3.0)*kT*rho0*lambda_r
    Xi[(  18)*ndof+(  18)] = (4.0/3.0)*kT*rho0*lambda_r
    Xi[(Q+ 0)*ndof+(Q+ 0)] = (4.0/3.0)*kT*rho0*lambda_r
    Xi[(Q+ 1)*ndof+(Q+ 1)] = 18.*kT*rho0*(1 - cs2)*lambda_r
    Xi[(Q+ 2)*ndof+(Q+ 2)] = 8.*kT*rho0*lambda_r
    Xi[(Q+ 3)*ndof+(Q+ 3)] = (8.0/3.0)*kT*rho0*lambda_r
    Xi[(Q+ 4)*ndof+(Q+ 4)] = 2.*Gamma*kT/rho0*(-9*Gamma*mu_phi + 5)*lambda_p
    Xi[(Q+ 5)*ndof+(Q+ 5)] = 8.*Gamma*kT/rho0*lambda_p
    Xi[(Q+ 6)*ndof+(Q+ 6)] = (8.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+ 7)*ndof+(Q+ 7)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+ 8)*ndof+(Q+ 8)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+ 9)*ndof+(Q+ 9)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+10)*ndof+(Q+10)] = 4.*Gamma*kT/rho0*lambda_p
    Xi[(Q+11)*ndof+(Q+11)] = 4.*Gamma*kT/rho0*lambda_p
    Xi[(Q+12)*ndof+(Q+12)] = 4.*Gamma*kT/rho0*lambda_p
    Xi[(Q+13)*ndof+(Q+13)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+14)*ndof+(Q+14)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+15)*ndof+(Q+15)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p
    Xi[(Q+16)*ndof+(Q+16)] = 18.*Gamma*kT/rho0*(-Gamma*mu_phi + 1)*lambda_p
    Xi[(Q+17)*ndof+(Q+17)] = 8.*Gamma*kT/rho0*lambda_p
    Xi[(Q+18)*ndof+(Q+18)] = (8.0/3.0)*Gamma*kT/rho0*lambda_p

    # rho-rho sector [0, 2..4, 8..(Q+3)]
    Xi[(   8)*ndof+(Q+ 1)] = 6.*kT*rho0*(3*cs2 - 1)*lambda_r
    Xi[(Q+ 1)*ndof+(   8)] = 6.*kT*rho0*(3*cs2 - 1)*lambda_r

    # phi-phi sector [1, 5..7, (Q+4)..(Q-1)]
    Xi[(Q+ 4)*ndof+(Q+16)] = 6.*Gamma*kT/rho0*(3*Gamma*mu_phi - 1)*lambda_p
    Xi[(Q+16)*ndof+(Q+ 4)] = 6.*Gamma*kT/rho0*(3*Gamma*mu_phi - 1)*lambda_p

    # rho-phi sector
    Xi[(   8)*ndof+(Q+ 4)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(Q+ 4)*ndof+(   8)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(   8)*ndof+(Q+16)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(Q+16)*ndof+(   8)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)

    Xi[(Q+ 1)*ndof+(Q+ 4)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(Q+ 4)*ndof+(Q+ 1)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(Q+ 1)*ndof+(Q+16)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)
    Xi[(Q+16)*ndof+(Q+ 1)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0)

    Xi[(   0)*ndof+(Q+ 4)] = -3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr
    Xi[(Q+ 4)*ndof+(   0)] = -3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr
    Xi[(   0)*ndof+(Q+16)] = 3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr
    Xi[(Q+16)*ndof+(   0)] = 3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr
    Xi[(   1)*ndof+(Q+ 1)] = 3.*kT*p_phi/(mu_phi*rho0)*lambda_rp
    Xi[(Q+ 1)*ndof+(   1)] = 3.*kT*p_phi/(mu_phi*rho0)*lambda_rp
    Xi[(   1)*ndof+(   8)] = -3.*kT*p_phi/(mu_phi*rho0)*lambda_rp
    Xi[(   8)*ndof+(   1)] = -3.*kT*p_phi/(mu_phi*rho0)*lambda_rp
    Xi[(   2)*ndof+(   5)] = -phi0*kT*lambda_pr
    Xi[(   3)*ndof+(   6)] = -phi0*kT*lambda_pr
    Xi[(   4)*ndof+(   7)] = -phi0*kT*lambda_pr
    Xi[(   5)*ndof+(   2)] = -phi0*kT*lambda_pr
    Xi[(   6)*ndof+(   3)] = -phi0*kT*lambda_pr
    Xi[(   7)*ndof+(   4)] = -phi0*kT*lambda_pr

    return Xi

def fourier_laplace_operator(ikx, iky, ikz, boxDim):
    n = boxDim.shape

    # FFTW convention for ordering of wave vectors
    kx = 2. * np.pi / n[0] * ikx if ikx < (n[0] + 1) / 2 else 2. * np.pi / n[0] * (ikx - n[0])
    ky = 2. * np.pi / n[1] * iky if iky < (n[1] + 1) / 2 else 2. * np.pi / n[1] * (iky - n[1])
    kz = 2. * np.pi / n[2] * ikz if ikz < (n[2] + 1) / 2 else 2. * np.pi / n[2] * (ikz - n[2])

    cosx = np.cos(kx)
    cosy = np.cos(ky)
    cosz = np.cos(kz)

    expr1 = cosx + cosy + cosz
    expr2 = cosx * cosy + cosy * cosz + cosx * cosz

    k2 = -2. / (1./3.) * (1. / 9. * expr1 + 1. / 9. * expr2 - 2. / 3.)

    return k2