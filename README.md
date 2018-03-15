# Mie Library

[![Build Status](https://travis-ci.org/darothen/py-mie.svg?branch=master)](https://travis-ci.org/darothen/py-mie) [![DOI](https://zenodo.org/badge/31327772.svg)](https://zenodo.org/badge/latestdoi/31327772)


This library wraps several Fortran implementations of Mie theory calculators ([1], [2], ... ) with a Python interface. The Fortran code itself has not been changed from their original form, and two-levels of drivers are available to run them. At the moment, an interface to the raw Mie codes (with all input parameters and output quantities) is not provided, although this functionality will probably be added in the future.

## Installation

Currently, you can install the package directly from the source code or via pip. To install directly from source, download the source code and execute:

``` bash
    $ pip install -e .
```

The package can also be installed directly from GitHub via the command:

``` bash
    $ pip install git+https://github.com/darothen/py-mie.git@branch-name
```


The only package dependency is ``f2py`` via ``numpy`` or ``scipy``.

## Unit Tests

To run the unit tests using nosetests, run the following:

``` bash
    $ nosetests
```


## Authors

This library was packaged by [Daniel Rothenberg (Massachusetts Institute of Technology)](http://www.github.com/darothen). The Mie code and specialized interfaces have been provided by:

- Rahul Zaveri (Pacific Northwest National Labs)
- Alexander Avramov (Massachusetts Institute of Technology)

## Citation

If this wrapper library was useful to you, please consider citing its [DOI](https://doi.org/10.5281/zenodo.192510) along with the relevant scientific citation to one of the **References** below. It's up to us to change our academic culture and ensure that researchers who take the time to build useful tools are properly credited for their hard work!

## References

1. [Owen B. Toon and T. P. Ackerman, "Algorithms for the calculation of scattering by stratified spheres," Appl. Opt. 20, 3657-3660 (1981)](http://dx.doi.org/10.1364/AO.20.003657)
2. [C. F. Bohren and D. R. Huffman, Absorption and Scattering of Light by Small Particles, New York, Wiley, 1983, 530 pages](http://onlinelibrary.wiley.com/book/10.1002/9783527618156)
