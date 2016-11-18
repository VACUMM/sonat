"""Compress and scale data

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
import numpy
import MV2
import cdms2
from vcmq import dict_filter, broadcast, get_atts, set_atts

from .__init__ import _Base_, PyARMError

npy = numpy
default_missing_value = npy.ma.default_fill_value(0.)


class Packer(_Base_):
    """Class to handle a single variable

    This class packs a variable in 2D space by removing
    all masked points and storing a space-time array.
    It performs this operation also on the weights.
    It is used for removing unnecessary points and
    simplifying the input format for analysis functions.

    Parameters
    ----------
    data: numpy or masked array
        An array whose first dimension is not compressible, like time
        or member.
    norm:
        Normalisation coefficients
    mask:
        Integer mask where valid data = 1
    """
    def __init__(self, data, norm=None, mean=None, nordim=None, logger=None, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Guess data type and copy
        self.input = input = data
        if cdms2.isVariable(data):
            self.array_type = 'MV2'
            self.array_mod = MV2
            if nordim is None:
                nordim = data.getTime() is None
            data = data.asma()
        elif npy.ma.isMA(data):
            self.array_type = 'numpy.ma'
            self.array_mod = numpy.ma
        else:
            self.array_type = 'numpy'
            self.array_mod = numpy
        self.data = data = data.copy().astype('d')
        self.dtype = data.dtype

        # Shape
        if nordim is None:
            self.logger.warning("Can't guess if data has a record axis "
                 "(like time). Please specify the record parameter. "
                 "Assumed without record axis.")
            nordim = True
        self.nordim = nordim
        self.shape = data.shape
        self.ndim = data.ndim
        if not self.nordim:
            self.nr = self.shape[0]
            self.nstot = data.size/self.nr
        else:
            self.nr = 0
            self.nstot = data.size
        self.nsdim = data.ndim - self.nrdim
        self.sshape = self.shape[self.nrdim:]

        # Store some info

        # - record
        if not self.nordim:
            if self.ismv2:
                self.raxis = input.shape[0]
            else:
                self.raxis = input.getAxis(0)
        else:
            self.raxis = None

        # - others axes and attributes
        if self.ismv2: # cdms -> ids

            self.saxes = input.getAxisList()[self.nrdim:]
            self.id = input.id
            self.atts =  {}
            for att in input.listattributes() + ['id']:
                val = getattr(input, att, att)
                self.atts[att] = val
                if att in ['units', 'long_name', 'id']:
                    setattr(self, att, val)
            self.grid = input.getGrid()

        else: # numpy/ma -> length

            self.saxes = data.shape[self.nrdim:]
            self.id = None
            self.atts = None
            self.grid = None

        # - missing value
        if self.isma:
            self.missing_value = data.get_fill_value()
        else:
            self.missing_value = default_missing_value
            data = numpy.ma.masked_values(data, default_missing_value)


        # Masking nans
        nans = npy.isnan(data)
        if nans.any():
            self.logger.warning("Masking %i NaNs"%nans.sum())
            if self.array_type == 'numpy':
                self.array_type = 'numpy.ma'
                self.array_mod = numpy.ma
                self.data = data = npy.ma.masked_where(nans, data, copy=False)
            else:
                data[nans] = npy.ma.masked

        # Mask
        self.good = ~npy.ma.getmaskarray(data)
        if self.withrdim:
            self.good = self.good.all(axis=0)
        self.ns = self.good.sum()
        self.compress = self.ns != self.good.size

        # Scale unpacked data
        self.masked = not self.good.any()
        if self.masked:
            self.logger.warning('all your data are masked')
            self.norm = 1.
            self.mean = 0
        else:
            # - mean
            if not self.withrdim:
                self.mean = 0
            elif mean is not None:
                self.mean = mean
            else:
                self.mean = data.mean(axis=0)
            # - normalisation factor
            if norm is True or norm is None:
                norm = data.std() # Standard norm
            elif norm is not False:
                if norm <0: # Relative norm, else strict norm
                    norm = abs(norm)*data.std()
            else:
                norm = 1.
            self.norm = norm
            # - apply
            self.scale(data)

        # Fill data
        data_num = data.filled(self.missing_value)

        # Pack
        self.packed_data = self.core_pack(data_num, force2d=False)

    @property
    def withrdim(self):
        return not self.nordim

    @property
    def nrdim(self):
        return int(not self.nordim)

    @property
    def ismv2(self):
        return self.array_type=='MV2'

    @property
    def isma(self):
        return  self.array_type=='numpy.ma'

    @property
    def isnumpy(self):
        return  self.array_type=='numpy'

    def core_pack(self, data_num, force2d=False):
        """Compress data along space if needed

        :Parameters:

        *data_num**: Pure numpy array.
        """
        # Check shape
        if data_num.shape[-self.nsdim:] != self.sshape:
            raise PyARMError("Data to pack has a wrong shape: {}".format(data_num.shape))

        # Remove bad channels ?
        nxdim = data_num.ndim-self.nsdim # dims other than space
        if self.compress: # Pack

            sl = [slice(None)]*nxdim+[self.good]
            pdata = data_num[sl].T

        else: # Don't pack, just reshape

            if nxdim+self.nsdim > 2:
                data_num.shape = data_num.shape[:nxdim]+(-1, )
            pdata = data_num.T

        # At least 1D
        pdata = npy.atleast_1d(pdata)

        # 2D?
        if force2d and pdata.ndim==1:
            if self.nsdim==0:
                pdata = npy.atleast_2d(pdata)
            else:
                pdata.shape = pdata.size, 1

        return npy.asfortranarray(pdata)

    def scale(self, data, copy=False, mean=None, norm=None, mode=None):
        """Remove mean and normalize unpacked data"""
        if mode=='mean':
            if norm is None:
                norm = False
            if mean is None:
                mean = True
        elif mode=='norm':
            if norm is None:
                norm = True
            if mean is None:
                mean = False
        if copy:
            data = data.clone() if cdms2.isVariable(data) else data.copy()
        if mean is not False:
            data[:] -= self.mean if mean is True or mean is None else mean
        if norm is not False:
            data[:] /= self.norm if norm is True or norm is None else norm
        return data

    def rescale(self, data, copy=False, mean=None, norm=None, mode=None):
        """Re-add mean and unnormalize unpacked data"""
        if mode=='mean':
            if norm is None:
                norm = False
            if mean is None:
                mean = True
        elif mode=='norm':
            if norm is None:
                norm = True
            if mean is None:
                mean = False
        if copy:
            data = data.clone() if cdms2.isVariable(data) else data.copy()
        if norm is not False:
            data[:] *= self.norm if norm is True or norm is None else norm
        if mean is not False:
            data[:] += self.mean if mean is True or mean is None else mean
        return data


    def repack(self, data, scale=True, force2d=False):
        """Pack a variable using previously computed mask array"""

        # Scale
        if scale:
            data = self.scale(data, copy=True)

        # Numpy
        if npy.ma.isMA(data):
            data = data.filled(default_missing_value)

        # Check shape
        nsdim = self.ndim-1
        if nsdim!=0:
            if data.shape[-nsdim:] != self.shape[-nsdim:]:
                self.error('Incompatible spatial shape (%s instead of %s)'
                    %(data.shape[-nsdim:], self.shape[-nsdim:]))
        # Pack it
        return self.core_pack(data, force2d=force2d)

    def _get_firstdims_(self, firstdims, firstaxes):
        shape = self.sshape

        if firstdims is not None:

            if not firstdims:
                firstdims = False
            elif not isinstance(firstdims, tuple):
                firstdims = (firstdims, )

        if firstdims is not False:

            if firstaxes is None and self.withrdim:
                firstaxes = [self.raxis]
            elif isinstance(firstaxes, tuple):
                firstaxes = list(firstaxes)
            elif not isinstance(firstaxes, list):
                firstaxes = [firstaxes]
            else:
                firstaxes = []

            if firstdims is None and firstaxes:
                firstdims = tuple([(isinstance(a, (int, long))  and a or len(a)) for a in firstaxes])
            shape = firstdims + shape
            if firstaxes and isinstance(firstaxes[0], int): # real axes, not ints
                firstaxes = None

        elif len(shape)==0:

            shape = (1, )

        return shape, firstdims, firstaxes

    def create_array(self, firstdims=None, format=1, firstaxes=None):
        """Initialize an array similar to input array

        Type of array is guessed from attribute :attr:`array_type`.

        The array is filled with attribute :attr:`missing_value` if
        pure numpy, else with masked values.

        :Params:

        *firstdims**, optional: Size of the first dimension(s).
              Defaults to attribute :attr:`nt`.
        *format**, optional: To format output array as a CDAT variable
              using information from analyzed array.

                - ``1``: Add all initial axes but the first one.
                - ``2`` or ``"full"``: Add attributes.
                - else, does not format.

        *firstaxes**, optional: Set the axes as the first ones.
              If ``firstaxes is None`` and ``firstdims is None``, it defaults
              to the first axis of analyzed array.
              You may also provide integers instead of axes, which are then
              used to set length of first axes if ``firstdims`` is not
              set.


        """
        # Get base array
        # - shape
        shape, firstdims, firstaxes = self._get_firstdims_(firstdims, firstaxes)
        # - create
        MM = eval(self.array_type)
        data = MM.zeros(shape, self.dtype)
        # - missing values
        if self.isnumpy:
            data[:] = npy.ma.masked
            data.set_fill_value(self.missing_value)
        else:
            data[:] = self.missing_value

        # Format CDAT variables
        if self.ismv2 and format:

            # Axes
            if not firstdims:
                firstdims = ()
            for i, axis in enumerate(self.saxes):
                data.setAxis(i+len(firstdims), axis)
            if firstdims is not False and firstaxes is not None:
                for i, a in enumerate(firstaxes):
                    if not isinstance(a, (int, long)):
                        try:
                            data.setAxis(i, a)
                        except:
                            pass

            # Attributes
            if format=='full' or format>1:
                set_atts(data, self.atts)
        return data


    def unpack(self, pdata, rescale=True, format=1, firstdims=None, firstaxes=None):
        """Unpack data along space, reshape, and optionally unnormalize and remean.

        Input is sub_space:other, output is other:split_space.

        Params
        ------
        pdata: 2D array
            Packed data.
        rescale: boolean, string
            Rescale (mean and norm) using :meth:`rescale`?
        format: int, "full"
            Format the variable using :meth:`create_array`?
        firstaxes: List of ints or 1d arrays
            First axes as used by :meth:`create_array`
        """


        # Unpack
        # - space is last
        pdata = npy.ascontiguousarray(pdata.T).copy()
        # - first dimensions
        if firstdims is None and firstaxes is None:
            firstdims = 0 if pdata.ndim==1 and not self.withrdim else pdata.shape[0]
        shape, firstdims, firstaxes = self._get_firstdims_(firstdims, firstaxes) #FIXME:remove _get_firstdims_?
        # - create variable
        data = self.create_array(firstdims=firstdims, format=format, firstaxes=firstaxes)
        # - check packed data shape
        firstdims = data.shape[:len(data.shape)-self.nsdim]
        if len(firstdims) > 1:
            pdata = pdata.reshape(firstdims + (-1, ))
        # - uncompress
        first_slices = (slice(None), )*len(firstdims) # max(1, len(firstdims)) fix for notime case
        if self.compress:
            mdata = data.asma() if self.ismv2 else data
            slices = first_slices+(self.good, )
            mdata[slices] = pdata # here?
            data[:] = mdata # just to be sure
        else:
            if self.nsdim==0 and pdata.ndim>len(firstdims):
                pdata = pdata[first_slices+(0, )]
            data[:] = pdata.reshape(data.shape)
        del pdata
        # - mask
        data[:] = npy.ma.masked_values(data, default_missing_value,  copy=0)

        # Rescale
        if rescale:
            self.rescale(data, mode=rescale)

        return data

    def get_record_axis(self, nr=None, offset=0):
        """Get the time axis or dimension of an input variable

        If CDAT is not used, length of axis is returned or nr if different.
        """
        if not self.withrdim:
            return
        axis = self.raxis
        if not isinstance(axis, int) and nr is not None and nr!=self.nr:
            dt = npy.median(npy.diff(axis[:]))
            axiso = cdms2.createAxis(npy.arange(nr)*dt+axis[0])
            for att, val in axis.attributes.items():
                setattr(axiso, att, val)
            axiso.id = axis.id
            axiso[:] += offset
            return axiso
        if not isinstance(axis, int) and offset:
            axis = axis.clone()
            axis[:] += offset
        elif isinstance(axis, int) and nr is not None and axis!=nr:
            return nr
        return axis

