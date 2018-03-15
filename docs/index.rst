.. py-mie documentation master file, created by
   sphinx-quickstart on Wed Mar 14 11:04:32 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to py-mie's documentation!
==================================

This library wraps several Fortran implementations of Mie theory calculators
([1], [2], ...) with a python interface. The Fortran code itself has not been
changed from their original form, and two-levels of drivers are available to run
them. At the moment, an interface to the raw Mie codes (with all input parameters
and output quantities) is not provided, although this functionality will probably
be added in the future.


Installation
------------

Currently, installation must be from the source code itself. Download the source
code and execute::

    $ python setup.py install

Alternatively, you can install directly from GitHub::

    $ pip install git+https://github.com/darothen/py-mie.git


Requirements
------------

The only package dependency is ``f2py`` via either ``numpy`` or ``scipy``.


Getting Help
-------------

To raise issues are ask specific questions about the library, please submit an
issue via the `GitHub repository <https://github.com/darothen/py-mie/issues>`_.


Running Unittests
-----------------

To run the unit tests using nosetests, run the following::

    $ nosetests


Authors
-------

This library was packaged by Daniel Rothenberg (Massachusetts Institute of Technology).
The Mie code and specialized interfaces have been provided by:

* Rahul Zaveri (Pacific Northwest National Labs)
* Alexander Avramov (Massachusetts Institute of Technology)


Citation
--------

If this wrapper library was useful to you, please consider citing its DOI along
with the relevant scientific citation to one of the References below. It's up to
us to change our academic culture and ensure that researchers who take the time
to build useful tools are properly credited for their hard work!

In addition, this package can be cited using the following DOI:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.192510.svg
   :target: https://doi.org/10.5281/zenodo.192510

----------
References
----------

1. Owen B. Toon and T.P. Ackerman, "Algorithms for the calculation of scattering by stratified spheres", Appl. Opt. 20, 3657-3660 (1981)
2. C.F. Bohren and D.R. Huffman, Absorption and Scattering of Light by Small Particles, New York, Wiley, 1983, 530 pages


=============
API Reference
=============

.. autofunction:: mie.bhmie_scatter
.. autofunction:: mie.core_shell_scatter
.. autofunction:: mie.integrate_mode
