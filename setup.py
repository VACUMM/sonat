import os
import sys
import re
import shutil
from numpy.distutils.core import setup, Extension
#from numpy.f2py import crackfortran

rootdir = os.path.dirname(__file__)

# Get values from the lib
keys = ['__version__', '__author__', '__date__']
f = open(os.path.join(rootdir, 'pyarm/__init__.py'))
for line in f:
    line = line[:-1].strip()
    for key in keys:
        if line.startswith(key):
            exec(line)
            break
f.close()

# Some info
version = __version__
date = __date__
description = 'Objective coastal observation array assement tools'
author = __author__
url = "http://www.ifremer.fr/pyarm"
license = "CeCILL"

# Special setups
if __name__=='__main__':

    # Config file
    cfg_file = os.path.join(rootdir, 'setup.cfg')
    default_cfg_file = os.path.join(rootdir, 'setup.default.cfg')
    if not os.path.exists(cfg_file):
        shutil.copy(default_cfg_file, cfg_file)


    # Setup the python module
    s = setup(name="pyarm",
        version = version,
        description = description,
        author = author,
        #author_email=author_email,
        maintainer = author,
        #maintainer_email=author_email,
        license = license,
        ext_modules = [
            Extension('pyarm._fcore', ['pyarm/fcore.f90'])
        ],
        #package_dir={'pyarm':'pyarm'},
        packages = ["pyarm"],

    )

