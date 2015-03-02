"""
Suite of tests to double check that the Fortran modules were built correctly
and are giving results in agreement with one another.
"""

from __future__ import absolute_import, division

import numpy as np
from numpy.testing import (assert_allclose, assert_array_less,
                           assert_equal, assert_raises)

from .._mie import ( bhmie_module, dmiess_module, dmilay_module )
#from _mie import bhmie_module, dmiess_module, dmilay_module

refr_BC = 1.95 + 0.79*1j
refr_OC = [
    real+imag*1j for real, imag in \
        zip([1.443, 1.420, 1.420, 1.420, 1.463, 1.510, 1.510, 1.520, 
             1.530, 1.530, 1.530, 1.530, 1.530, 1.124], 
            [0.006, 0.018, 0.011, 0.008, 0.016, 0.022, 0.019, 0.016, 
             0.007, 0.006, 0.005, 0.008, 0.030, 0.079])
]
refr_SO4 = [
    1.53+imag*1j for imag in \
        [0.158, 0.057, 0.003, 0.001, 0.001, 0.000, 0.000, 0.000, 
         0.000, 0.000, 0.000, 0.000, 0.000, 0.551]            
]

wv_max = [3.846, 3.077, 2.500, 2.151, 1.942, 1.626, 1.299, 1.242, 
          0.778, 0.625, 0.442, 0.345, 0.263, 12.195]
wv_min = [3.077, 2.500, 2.151, 1.942, 1.626, 1.299, 1.242,  .778, 
          .625, .442, .345, .263, .200,  3.846]
wavelengths = 0.5*(np.array(wv_max)+np.array(wv_min))


radii = [0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 
         0.25, 0.5, 1., 2.5, 5., 10., 25., 50., ]

core_fractions = np.linspace(0, 1, 11)[1:-1]

species_refr_inds = {
    'SO4': refr_SO4, 'BC': refr_BC, 'OC': refr_OC, 
}


def test_homogeneous(rtol=1.0):

    ERR_MSG = """
    {0:s} vs {1:s}
    ----------------------------------------------------
        wavelength: {2:>4.2f}
            radius: {3:>4.2f}
          refr_ind: {4:}
    """

    def _run_homogeneous(radius, refr_ind, wavelength):
        return [ bhmie_module.bhmie_driver(radius, refr_BC, wavelength),
                 dmiess_module.dmiess_driver(radius, radius, refr_BC, refr_BC, wavelength),
                 dmilay_module.dmilay_driver(radius, radius, refr_BC, refr_BC, wavelength) ]

    def _check_mie(species, radius, refr_ind, wavelength, rtol=1.0):

        err_msg_run = lambda l, r : ERR_MSG.format(l, r, wavelength, radius, refr_ind)

        mie_results = _run_homogeneous(radius, refr_ind, wavelength)
        bhmie  = np.array(mie_results)[0, :-1]
        dmiess = np.array(mie_results)[1, :-1]
        dmilay = np.array(mie_results)[2, :-1]

        assert_allclose(bhmie, dmiess, rtol=rtol, 
                        err_msg=err_msg_run('BHMIE', 'DMIESS'))
        assert_allclose(bhmie, dmilay, rtol=rtol, 
                        err_msg=err_msg_run('BHMIE', 'DMILAY'))
        assert_allclose(dmilay, dmiess, rtol=rtol, 
                        err_msg=err_msg_run('DMILAY', 'DMIESS'))

    for species in species_refr_inds:
        refr_inds = species_refr_inds[species]

        if species == "BC":
            for wavelength in wavelengths:
                for radius in radii:
                    yield _check_mie, species, radius, refr_BC, wavelength, rtol
        else:
            for (refr_ind, wavelength) in zip(refr_inds, wavelengths):
                for radius in radii:
                    yield _check_mie, species, radius, refr_ind, wavelength, rtol

def test_heterogeneous(rtol=1.0):
    """ Black carbon core, sulfate shell - DMIESS vs DMILAY"""

    ERR_MSG = """
    {0:s} vs {1:s}
    ----------------------------------------------------
        wavelength: {2:>4.2f}
            radius: {3:>4.2f}
         core_frac: {4:>4.1f}
    """

    def _run_heterogeneous(radius_shell, radius_core,
                           refr_shell, refr_core, wavelength):
        return [ dmiess_module.dmiess_driver(radius_shell, radius_core, 
                                             refr_shell, refr_core, wavelength),
                 dmilay_module.dmilay_driver(radius_shell, radius_core, 
                                             refr_shell, refr_core, wavelength) ]

    def _check_mie(radius, core_fraction, refr_ind_index, 
                   wavelength, rtol=1.0):

        radius_shell = radius
        radius_core = radius*core_fraction

        refr_shell = refr_SO4[refr_ind_index]
        refr_core = refr_BC

        err_msg_run = lambda l, r : ERR_MSG.format(l, r, wavelength, radius,
                                                   core_fraction)

        mie_results = _run_heterogeneous(radius_shell, radius_core,
                                         refr_shell, refr_core, wavelength)
        dmiess = np.array(mie_results)[0, :-1]
        dmilay = np.array(mie_results)[1, :-1]

        assert_allclose(dmilay, dmiess, rtol=rtol, 
                        err_msg=err_msg_run('DMILAY', 'DMIESS'))

    for (i, wavelength) in enumerate(wavelengths):
        for radius in radii:
            for core_fraction in core_fractions:
                yield _check_mie, radius, core_fraction, i, wavelength, rtol
