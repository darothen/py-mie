#!/usr/bin/env python

import os
from numpy.distutils.core import setup, Extension

src_dir = "src"
pkg_dir = "mie"

DESCRIPTION = "Python wrapper for Mie theory calculators"
NAME = "py-mie"
AUTHOR = "Daniel Rothenberg"
AUTHOR_EMAIL = "darothen@mit.edu"
MAINTAINER = "Daniel Rothenberg"
MAINTAINER_EMAIL = "darothen@mit.edu"
URL = 'http://github.com/darothen/py-mie'
DOWNLOAD_URL = 'http://github.com/darothen/py-mie'
LICENSE = 'MIT'

#import mie
#VERSION = mie.__version__

## Setup the compiled library modules
mie_files = ["mod_%s.f90" % s for s in ["kinds", "dmiess", "bhmie", "dmilay"]]
sources   = [os.path.join(src_dir, s) for s in mie_files] 
sources  += [os.path.join(pkg_dir, "mod_mie.pyf"), ] 

mie_ext = Extension(
    'mie._mie', sources=sources
)

setup(
    name=NAME,
    version='0.3.0',
    description=DESCRIPTION,
    #long_description=LONG_DESCRIPTION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    url=URL,
    download_url=DOWNLOAD_URL,
    license=LICENSE,
    #packages=['supersmoother',
    #          'supersmoother.tests',
    #      ],
    install_requires=['numpy', ],
    packages=['mie','mie.tests', ],
    classifiers=[
        'Development Status :: 3 - Alpha', 
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],

    ext_modules = [
        mie_ext, 
    ]
)
