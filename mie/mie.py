import numpy as np

from ._mie import dmiess_module as dmiess
from ._mie import bhmie_module as bhmie
from ._mie import dmilay_module as dmilay

#: Density of pure water, in kg/m^3
RHO_WATER = 1e3

__all__ = ["bhmie_scatter", "core_shell_scatter", "integrate_mode", ]


def bhmie_scatter(particle_radius, radiation_lambda, n_particle):
    """ Compute the scattering/absorption efficiency and asymmetry
    parameter for a homogenous particle.

    This function interfaces with the compiled Mie theory modules in order
    to determine the scattering parameters.

    Parameters
    ----------
    particle_radius : float
        The total particle radius (core + shell) in microns
    radiation_lambda : float
        Wavelength of incident radiation in microns
    n_particle : complex
        Complex refractive indices of the particle material

    Returns
    -------
    Qsca, Qabs, asym : floats
        The scattering efficiency, absorption efficiency, and asymmetry
        parameter for the specified particle

    Examples
    --------

    Evaluate the scattering for a PSL with particle radius of 0.5 microns and RI=1.5
    at an incident wavelength of 658 nm:

    >>> qsca, _, _ = mie.bhmie_scatter(0.5, 0.658, 1.5+1j*0)

    """

    # Pass directly to Mie module
    Qext0, Qsca0, asym0 = bhmie.bhmie_driver(particle_radius, n_particle,
                                             radiation_lambda)

    # Post-process to properly set scattering and absorption efficiencies
    Qsca = np.min([Qsca0, Qext0])  # scattering efficiency
    Qabs = Qext0 - Qsca            # absorption efficiency
    asym = asym0

    return Qsca, Qabs, asym


def core_shell_scatter(particle_radius, core_fraction, radiation_lambda,
                       n_shell, n_core):
    """ Compute the scattering/absorption efficiency and asymmetry
    parameter for a heterogeneous, core-shell mixed particle.

    This function interfaces with the compiled Mie theory modules in order
    to determine the scattering parameters.

    Parameters
    ----------
    particle_radius : float
        The total particle radius (core + shell) in microns
    core_fraction : float
        The fraction of the particle comprised by its core, 0.0-1.0
    radiation_lambda : float
        Wavelength of incident radiation in microns
    n_shell, n_core : complex
        Complex refractive indices of the shell, and core respectively

    Returns
    -------
    Qsca, Qabs, asym : floats
        The scattering efficiency, absorption efficiency, and asymmetry
        parameter for the specified particle

    """

    # Pass directly to Mie module
    Qext0, Qsca0, asym0 = dmiess.dmiess_driver(
        particle_radius,
        core_fraction*particle_radius,
        n_shell,
        n_core,
        radiation_lambda)

    # Post-process to properly set scattering and absorption efficiencies
    Qsca = np.min([Qsca0, Qext0])  # scattering efficiency
    Qabs = Qext0 - Qsca            # absorption efficiency
    asym = asym0

    return Qsca, Qabs, asym


def integrate_mode(core_fraction, n_shell, n_core, radiation_lambda,
                   mode_radius, mode_sigma,
                   r_min=1e-3, r_max=100., nr=200):
    """ Integrate Mie theory calculations over a lognormal aerosol mode with
    homogeneous particle properties, weighting by size distribution.

    Parameters
    ----------
    core_fraction : float
        The fraction of the particle comprised by its core, 0.0-1.0
    n_shell, n_core : complex
        Complex refractive indices of the shell, and core respectively
    radiation_lambda : float
        Wavelength of incident radiation in microns
    mode_radius : float
        The geometric mean or mode radius of the aerosol size distribution, in
        microns
    mode_sigma : float
        The geometric standard deviation of the aerosol size distribution
    r_min, r_max : float (optional)
        The minimum and maximum particle radii to use in the integration, in
        microns
    nr : int (optional)
        The number of particle radii to use in the integration

    Returns
    -------
    mie_sca, mie_abs, mie_asym : floats
        Scattering efficiency, absorption efficiency, and asymmetry parameter
        integrated over an aerosol size distribution

    """

    # Generate the integration grid for the particle size distribution
    logr, dlogr = np.linspace(np.log(r_min), np.log(r_max), nr, retstep=True)
    radii = np.exp(logr)

    sumsca  = 0.0
    sumabs  = 0.0
    sumg    = 0.0
    volwet  = 0.0
    volcore = 0.0

    # Integration loop
    for i, radius in enumerate(radii):

        # Mie theory calculation
        Qsca, Qabs, asym = core_shell_scatter(
            radius, core_fraction, radiation_lambda, n_shell, n_core
        )

        # Compute weights and volumes for integral sum
        exparg   = np.log(radius / mode_radius)/np.log(mode_sigma)
        dsdlogr  = np.exp(-0.5*exparg**2)  # m^2/m3(air)] log-normal cross section area
        volwet  += (4./3.)*(radius*1e-6)*dsdlogr*dlogr  # [m3/m3(air)] wet volume
        volcore += (4./3.)*(core_fraction**3)*(radius*1e-6)*dsdlogr*dlogr

        sumabs += Qabs*dsdlogr*dlogr  # [m^2/m3(air)] absorption cross-section
        sumsca += Qsca*dsdlogr*dlogr  # [m^2/m3(air)] scattering cross-section
        sumg   += asym*Qsca*dsdlogr*dlogr

    mie_sca  = np.max([sumsca, 0.]) / (volwet * RHO_WATER)  # [m^2/kg] specific scattering cross-section
    mie_abs  = np.max([sumabs, 0.]) / (volwet * RHO_WATER)  # [m^2/kg] specific absorption cross-section
    mie_asym =                 sumg / sumsca                # average asymmetry

    return mie_sca, mie_abs, mie_asym
