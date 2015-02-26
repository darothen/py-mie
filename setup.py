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
#LICENSE = 'BSD 3-clause'

#import mie
#VERSION = mie.__version__

## Setup the compiled library modules
sources = []

mie_ext = Extension(
    'mie._mie', sources=[os.path.join(src_dir, "mod_core_shell.f90"), 
                         os.path.join(pkg_dir, "mod_core_shell.pyf"), ]
)

setup(
    name=NAME,
    #version=VERSION,
    description=DESCRIPTION,
    #long_description=LONG_DESCRIPTION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    url=URL,
    download_url=DOWNLOAD_URL,
    #license=LICENSE,
    #packages=['supersmoother',
    #          'supersmoother.tests',
    #      ],
    packages=['mie',],
    classifiers=[
        #'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        #'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Fortran :: 90',
    ],

    ext_modules = [
        mie_ext, 
    ]
)
