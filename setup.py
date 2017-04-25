#!/usr/bin/env python
#
# Copyright IFREMER (2016-2017)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
import os
import sys
import re
import shutil
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

rootdir = os.path.dirname(__file__)

# Get values from the sources
keys = ['__version__', '__author__', '__date__',  '__url__', '__email__']
f = open(os.path.join(rootdir, 'sonat/__init__.py'))
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
description = 'Stochastic ocean Observing Network Assessment Toolkit'
author = __author__
author_email = __email__
url = __url__
license = "CeCILL"
long_description = open('README.rst').read()

# Special setups
if __name__=='__main__':

    # Config file
    cfg_file = os.path.join(rootdir, 'setup.cfg')
    default_cfg_file = os.path.join(rootdir, 'setup.default.cfg')
    if not os.path.exists(cfg_file):
        shutil.copy(default_cfg_file, cfg_file)

    # Distutils config
    def configuration(parent_package='',top_path=None):
        config = Configuration()
        config.add_scripts('bin/sonat') # main executable
        config.add_data_dir(('sonat/data', 'data')) # all data files
        return config

    # Setup the python module
    s = setup(name="sonat",
        version = version,
        description = description,
        long_description = long_description,
        author = author,
        author_email=author_email,
        maintainer = author,
        maintainer_email=author_email,
        license = license,
        ext_modules = [
            Extension('sonat._fcore', ['sonat/fcore.f90'])
        ],
        packages = ["sonat", "sonat.test"],
        package_dir = {"sonat.test": "test"},
        package_data = {"sonat": ["matplotlibrc", "sonat.ini",
                                  "CollapsibleLists.compressed.js",
                                  "sonat.css",
                                  "button-*.png", "list-item*.png"],
            "sonat.test": ["sonat.cfg", "sangoma*", "inputs/*.txt"]},
        configuration = configuration,

    )

