Mie Library
===========

This library wraps several implementations of Mie theory calculators ([1] ... )
with a Python interface. 

Installation
------------

Currently, installation must be from the source code itself. Download the source code and execute:

    $ python setup.py install

The only package dependency is ``f2py`` via ``numpy`` or ``scipy``.

Authors
-------

This library was developed by [Daniel Rothenberg][http://www.github.com/darothen]. The Mie code and specialized interfaces have been provided by:

- Rahul Zaveri (Pacific Northwest National Labs)
- Alexander Avramov (Massachusetts Institute of Technology)

References
----------

1. [Owen B. Toon and T. P. Ackerman, "Algorithms for the calculation of scattering by stratified spheres," Appl. Opt. 20, 3657-3660 (1981)](http://dx.doi.org/10.1364/AO.20.003657)
