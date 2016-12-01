"""Stacking module

Inspired from the spanlib library (http://www.github.com/stefraynaud/spanlib)
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

import numpy as npy
from vcmq import broadcast

from .__init__ import SONATError
from .misc import _Base_, _MapIO_
from .pack import Packer, default_missing_value

class _MapIO_(object):

    def _load_input_(self, input):
        self.input = input
        if isinstance(input, (list, tuple)):
            inputs = list(input)
            self.map = len(input)
        else:
            inputs = [input]
            self.map = 0
        return inputs

    def unmap(self, out):
        """Remap out so as to be in the same form as input data

        It has the opposite effect of :meth:`remap`.
        """
        if self.map==0:
            return out[0]
        return out

    def remap(self, values, reshape=True, fill_value=None):
        """Makes sure that values is a list (or tuple) of length :attr:`nvar`

        :func:`broadcast` is called when reshaping.

        It has the opposite effect of :meth:`unmap`.
        """
        # We always need a list
        if not isinstance(values, (list, tuple)):
            values = [values]

        # Check length
        n = len(self)
        if len(values)!=n:
            if not reshape:
                raise SONATError('Wrong number of input items (%i instead of %i)'
                    %(len(values), n))
            values = broadcast(values, n, fill_value)
        return values


class Stacker(_Base_, _MapIO_):
    """Class to handle one variable or a list of variables

    This fonction concatenates several dataset that have the
    same record axis. It is useful for analysing for example
    several variables at the same record.
    It takes into account weights, masks and axes.

    Parameters
    ----------
    input: single or list of arrays
        They must all have the same record length if any.
    """

    def __init__(self, input, norms=None, means=None, nordim=None, logger=None, **kwargs):

        # Logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Input
        self.inputs = inputs = _MapIO_._load_input_(input)
        self.nvar = self.nd = len(inputs)

        # Other inits
        self.packers = []
        self.datas = []
        self.nr = None
        norms = self.remap(norms, reshape=True)
        means = self.remap(means, reshape=True)
        if self.nvar==1 and norms[0] is None: norms = [False]
        self.masked = False
        self.nordim = nordim

        # Loop on inputs
        for idata, data in enumerate(inputs):

            # Create the Packer instance and pack array
            packer = Packer(data, norm=norms[idata], mean=means[idata], nordim=nordim)
            self.packers.append(packer)
            self.masked |= packer.masked
            self.datas.append(packer.data)

            # Check nr
            if self.nr is None:
                self.nr = packer.nr
            elif self.nr != packer.nr:
                raise SONATError(('Record dimension of variable {idata} must '
                    'be {self.nr} (not {packer.nr})').format(**locals()))

        # Merge
        self.stacked_data = npy.asfortranarray(
            npy.concatenate([p.packed_data for p in self.packers], axis=0))
        self.splits = npy.cumsum([p.packed_data.shape[0] for p in self.packers[:-1]])
        self.ns = self.stacked_data.shape[0]
        self.nrv = (self.stacked_data!=default_missing_value).any(axis=0).sum()

    def __len__(self):
        return self.nvar

    def __getitem__(self, i):
        return self.packers[i]

    @property
    def data(self):
        return self.unmap(self.datas)

    def has_cdat(self, idata=None):
        """Check if there were at least one input array of CDAT type (:mod:`MV2`)

        :Sea also: :meth:`Data.has_cdat`
        """
        # Single var
        if idata is not None:
            return self[idata].has_cdat()

        # At least one
        for d in self.packers:
            if d.has_cdat(): return True
        return False

    def get_record_axis(self, nr=None, offset=0, idata=0):
        """Get the record axis of an input variable

        If CDAT is not used, length of axis is returned.
        """
        return self[idata].get_record_axis(nr, offset)


    def get_norms(self, idata=None):
        """Get :attr:`norms` for one or all input variables"""
        if idata is None:
            return [p.norm for p in self.packers]
        return self[idata].norm
    norms = property(get_norms, doc="Normalization factors")

    def get_means(self, idata=None):
        """Get :attr:`means` for one or all input variables"""
        if idata is None:
            return [p.mean for p in self.packers]
        return self[idata].mean
    means = property(get_means, doc="Record averages")

    @property
    def norm(self):
        self.unmap(self.norms)

    @property
    def mean(self):
        self.unmap(self.means)

    def restack(self, input, scale=True):
        """Stack new variables as a fortran array

        It has the opposite effect of :meth:`unstack`.

        :Params:

            - **input**: Argument in the same form as initialization data.
            - **scale**, optional: Scale the variable (mean and norm), and optionally
              remove spatial mean if == 2.

        :Seel also: :meth:`Pack.repack`
        """

        # Check input data
        inputs = self.remap(input, reshape=False)

        # Pack
        packs = [self[idata].repack(data, scale=scale)
            for idata, data in enumerate(inputs)]

        # Check record length (first axis)
        nr1 = inputs[0].rsize / self[0].nstot
        if len(inputs)>1:
            for i, d in enumerate(inputs[1:]):
                i += 1
                nr2 = d.size/self[i].nstot
                if packs[0].ndim != packs[i].ndim:
                    if packs[0].ndim==1:
                        SONATError('Variable {} must not have a record dimension '.format(i))
                    else:
                        SONATError('Variable {} must have a record dimension '.format(i))
                elif nr1!=nr2:
                    raise SONATError(('Record length of variable {i} ({nr2}) '
                        'different from that of first variable ({nr1})').format(**locals()))

        # Stack
        sdata = npy.asfortranarray(npy.concatenate(packs, axis=0))

        return sdata


    def unstack(self, sdata, rescale=True, format=1, firstdims=None, **kwargs):
        """Unstack and unpack data

        It has the opposite effect of :meth:`restack`.


        :Params:

            - **sdata**: Stacked data (as computed by :meth:`__init__` or :meth:`restack`).
            - **rescale**, optional: Rescale the variable (mean and norm) and
              add spatial mean if == 2.
            - **format**, optional: Format the variable (see :meth:`Data.create_array`).
            - **firsdims**, optional: First axis (see :meth:`Data.create_array`).


        :Seel also: :meth:`Data.unpack`
        """

        # Unstack
        packs = npy.split(sdata, self.splits, axis=0)

        # Unpack
        return [self[i].unpack(pdata, rescale=rescale, format=format,
            firstdims=firstdims, **kwargs) for i, pdata in enumerate(packs)]


    def create_arrays(self, **kwargs):
        """Create arrays similar to input arrays"""
        return [p.create_array(**kwargs) for p in self.packers]

    def format_arrays(self, data, **kwargs):
        """Format an array as input"""
        data = self.remap(data)
        return [self[i].format_array(d, **kwargs) for i, d in enumerate(data)]

    def fill_arrays(self, data, **kwargs):
        """Fill arrays similar to input arrays with data"""
        data = self.remap(data)
        arrs = self.create_arrays(**kwargs)
        for i, d in enumerate(data):
            arrs[i][:] = d
        return arrs
