"""
User defined stuff.

See :func:`load_user_code_file` for more information.
"""
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
import imp
import inspect

from vcmq import Dataset, register_dataset #, register_cf_variable

from .__init__ import SONATError, sonat_warn
from .obs import register_obs_platform, NcObsPlatform
from .arm import (register_arm_score_function, ARM_SCORE_FUNCTION_PREFIX,
                  ARMSA, register_arm_sensitivity_analyser)

#: Default user defined code file name
SONAT_USER_CODE_FILE = 'mysonat.py'


def load_user_code_file(myfile=None):
    """Load a user code file as a module name :mod:`mysonat` and
    register its stuff

    Observation platforms
        It registers as new observation platform all subclasses of
        :class:`sonat.obs.NcObsPlatform` with a valid :
        :attr:`platform_type` attribute, using function
        :func:`~sonat.obs.register_platform`.
    ARM score functions
        It registers as new ARM score functions all functions whose name
        starts with :attr:`~sonat.arm.ARM_SCORE_FUNCTION_PREFIX` and
        who takes 3 mandatory arguments (ev, arm, rep), using the
        :func:`~sonat.arm.register_arm_score_function`.
    ARM sensitivity analysers
        It registers as new ARM sensitivity analysers all classes that
        are subclasses of :attr:`~sonat.arm.ARMSA`.
    CF variables
        It registers new standard :mod:`vacumm.data.cf` variables
        declared in the dictionary :attr:`vacumm_cf_variables`, with
        the key as the short name and the value as specifications passed
        as arguments to :func:`vacumm.data.cf.register_cf_variable`.
    Generic datasets
        It registers all declared subclassed of
        :class:`vacumm.data.misc.dataset.Dataset` using the
        :func:`vacumm.data.register_dataset` function.
    """

    # Load module
    if myfile is None:
        myfile = SONAT_USER_CODE_FILE
        warn = False
    else:
        warn = True
    if not os.path.exists(myfile):
        if warn:
            sonat_warn('User code file not found: ' + myfile)
        return
    mysonat = imp.load_source('mysonat', myfile)

    # Loop on members
    for member in inspect.getmembers(mysonat):

        # Local only
        if inspect.getmodule(member[1]) is not mysonat:
            continue
        name, obj = member

        # User defined platforms
        if issubclass(obj, NcObsPlatform):

            register_obs_platform(obj)

        # User defined arm sensitivity analyser
        if issubclass(obj, ARMSA):

            register_arm_sensitivity_analyser(obj)

        # User defined datasets
        elif issubclass(obj, Dataset):

            register_dataset(obj)

        # User defined vacumm.cf variables
        elif False and name.lower() == 'vacumm_cf_variables':

            # Check type
            if not isinstance(obj, dict):
                raise SONATError('vacumm_cf_variables must be dictionary')

            # Loop on specs
            for varname, specs in obj.items():

                register_cf_variable(varname, **specs)

        # User defined ARM score functions
        elif name.lower().startswith(ARM_SCORE_FUNCTION_PREFIX):

            register_arm_score_function(name, obj)
