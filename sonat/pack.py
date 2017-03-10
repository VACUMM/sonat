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
from vcmq import dict_filter, broadcast, get_atts, set_atts, isaxis

from .__init__ import SONATError
from .misc import _Base_

N = npy = numpy
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
    def __init__(self, data, norm=None, mean=None, nordim=None, logger=None,
            missing_value=None, **kwargs):

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
            self.rshape = input.shape[:1]
            self.rsize = N.multiply.reduce(self.rshape)
        else:
            self.raxis = None
            self.rshape = ()
            self.rsize = 0
        self.nr = self.rsize

        # - others axes
        if self.ismv2: # cdms -> ids

            self.saxes = input.getAxisList()[self.nrdim:]
            self.grid = input.getGrid()

        else: # numpy/ma -> length

            self.saxes = data.shape[self.nrdim:]
            self.grid = None

        # - missing value
        if self.isma:
            if missing_value is not None:
                data.set_fill_value(missing_value)
            self._missing_value = data.get_fill_value()
        else:
            if missing_value is None:
                missing_value = default_missing_value
            self._missing_value = missing_value
            data = numpy.ma.masked_values(data, self._missing_value)


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

        # Scale unpacked data
        if self.masked:
            self.logger.warning('all your data are masked')
            self._norm = 1.
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
                norm = (data-self.mean).std() # Standard norm
            elif norm is not False:
                if norm <0: # Relative norm, else strict norm
                    norm = abs(norm)*(data-self.mean).std()
            else:
                norm = 1.
            self._norm = norm
            # - apply
            self.scale(data)

        # Fill data
        data_num = data.filled(self._missing_value)

        # Pack
        self.packed_data = self.core_pack(data_num, force2d=False)

    @property
    def ns(self):
        return self.good.sum()

    @property
    def compress(self):
        return self.ns != self.good.size

    @property
    def masked(self):
        return not self.good.any()

    @property
    def psize(self):
        return self.packed_data.size

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

    @property
    def id(self):
        if self.ismv2:
            return self.input.id

    varname = id

    @property
    def atts(self):
        if not self.ismv2:
            return
        atts = {}
        for att in self.input.listattributes() + ['id']:
            val = getattr(self.input, att)
            atts[att] = val
        return atts

    @property
    def units(self):
        if not self.ismv2:
            return
        return self.atts.get('units', None)

    @property
    def long_name(self):
        if not self.ismv2:
            return
        return self.atts.get('long_name', None)

    def set_norm(self, norm):
        self.packed_data *=  self._norm / norm
        self._norm = norm

    def get_norm(self):
        return self._norm

    norm = property(fget=get_norm, fset=set_norm, doc="Normalisation factor")


    def set_missing_value(self, missing_value):
        if missing_value != self._missing_value:
            self.packed_data[N.isclose(self.packed_data, self._missing_value)] = \
                missing_value
            self._missing_value = missing_value

    def get_missing_value(self):
        return self._missing_value

    missing_value = property(fget=get_missing_value, fset=set_missing_value)
    fill_value = missing_value


    def core_pack(self, data_num, force2d=False):
        """Compress data along space if needed

        :Parameters:

        *data_num**: Pure numpy array.
        """
        # Check shape
        if data_num.shape[-self.nsdim:] != self.sshape:
            raise SONATError("Data to pack has a wrong spatial shape: {} (!= {})".format(
                data_num.shape[-self.nsdim:], self.sshape))

        # Remove bad channels ?
        nxdim = data_num.ndim-self.nsdim # dims other than space
        if self.compress: # Pack

            sl = [slice(None)]*nxdim+[self.good]
            pdata = data_num[sl].T

        else: # Don't pack, just reshape

            if self.nsdim > 1:
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
            data[:] /= self._norm if norm is True or norm is None else norm
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
        elif mode is False:
            mean = norm = False
        if copy:
            data = data.clone() if cdms2.isVariable(data) else data.copy()
        if norm is not False:
            data[:] *= self._norm if norm is True or norm is None else norm
        if mean is not False:
            data[:] += self.mean if mean is True or mean is None else mean
        return data


    def repack(self, data, scale=True, force2d=False, missing_value=None):
        """Pack a variable using previously computed mask array"""

        # Scale
        if scale:
            data = self.scale(data, copy=True)

        # Numpy
        if npy.ma.isMA(data):
            if missing_value is None:
                missing_value = self._missing_value
            data = data.filled(missing_value)

        # Check shape
        nsdim = self.ndim-1
        if nsdim!=0:
            if data.shape[-nsdim:] != self.shape[-nsdim:]:
                self.error('Incompatible spatial shape (%s instead of %s)'
                    %(data.shape[-nsdim:], self.shape[-nsdim:]))
        # Pack it
        return self.core_pack(data, force2d=force2d)

    def _get_firstdims_(self, firstdims, pshape=None):
        """Get consistent first dims and axes

        Parameters
        ----------
        firstdims: int, tuple, None, True, False
            Dimension sizes or axes
        pshape: None, tuple
            Packed shape with spatial dimension in last position

        """

        # With record dim?
        if pshape is not None:
            if len(pshape)==1:
                firstdims = False

        # No first dims, only space
        if firstdims is False: # no record dim
            return (), []

        # From input
        if firstdims is None or firstdims is True:

            # We really want it as the input so we check
            if pshape is not None and firstdims is True and self.rsize!=pshape[0]:
                raise SONATError("Wrong requested shape")

            # Ok, no chek or valid -> like input
            if pshape is None or self.rsize==pshape[0]:
                return self.rshape, [self.raxis] if self.withrdim else []

        # Explicit first dims?
        if firstdims is not None:

            # Transform to pure axes and dims
            if not isinstance(firstdims, tuple):
                firstdims = (firstdims, )
            firstdims_ = ()
            firstaxes = []
            for dim in firstdims:
                if isinstance(dim, int):
                    firstdims_ += dim,
                    firstaxes.append(None)
                else:
                    firstdims_ += len(dim),
                    firstaxes.append(dim)

            # Check
            if pshape is not None and N.multiply.reduce(firstdims_)!=pshape[0]:
                raise SONATError("Wrong requested shape")

            return firstdims_, firstaxes

        # From pshape only
        return pshape[:1], [None]

    def create_array(self, firstdims=None, format=1, pshape=None, id=None,
                     atts=None, format_atts=True):
        """Initialize an array similar to input array

        Type of array is guessed from attribute :attr:`array_type`.

        The array is filled with attribute :attr:`missing_value` if
        pure numpy, else with masked values.

        Parameters
        ----------
        firstdims: optional
            Size of the first dimension(s) or list of axes.
            Defaults to attribute :attr:`nt`.
        format: optional
            To format output array as a CDAT variable
            using information from analyzed array.

                - ``1``: Add all initial axes but the first one.
                - ``2`` or ``"full"``: Add attributes.
                - else, does not format.


        """
        # Get base array
        # - shape
        firstdims, firstaxes = self._get_firstdims_(firstdims, pshape=pshape)
        # - create
        MM = eval(self.array_type)
        data = MM.zeros(firstdims + self.sshape, self.dtype)
        # - missing values
        if self.isnumpy:
            data[:] = npy.ma.masked
            data.set_fill_value(self._missing_value)
        else:
            data[:] = self._missing_value

        # Format CDAT variables
        if self.ismv2 and format:

            data = self._format_array_(data, firstdims, firstaxes, format,
                id=id, atts=atts, format_atts=format_atts)

        return data

    def format_array(self, data, mode=1, firstdims=None, id=None, atts=None,
                     format_atts=True):
        """Format to a MV2 array similar to input array"""
        # Input was not formatted
        if not self.ismv2:
            return data

        # Make sure to have a MV2 array
        data = MV2.asarray(data)

        # Record or other first dims
        rdims = data.shape[:-self.nsdim]
        pshape = (N.multiply.reduce(rdims), ) if rdims else ()
        pshape += 1,  # fake spatial dim

        # First dimes
        firstdims, firstaxes = self._get_firstdims_(firstdims, pshape=pshape)

        # Format
        return self._format_array_(data, firstdims, firstaxes, mode,
            id=id, atts=atts, format_atts=format_atts)

    def _format_array_(self, data, firstdims, firstaxes, mode, id=None, atts=None,
                       format_atts=True):

        if firstaxes:
            for i, a in enumerate(firstaxes):
                if a is not None and not isinstance(a, int):
                    try:
                        data.setAxis(i, a)
                    except:
                        pass

        # Space
        for i, axis in enumerate(self.saxes):
            data.setAxis(i+len(firstdims), axis)

        # Id
        if atts is None:
            atts = self.atts.copy()
            atts['id'] = self.id
        else:
            mode = 2
        if id is True:
            id = self.id
        elif isinstance(id, basestring):
            id = id.format(**self.atts)
        else:
            id = None
        if id:
            data.id = id

        # Attributes
        if mode=='full' or mode>1:
            if 'id' in atts:
                del atts['id']

            # Format string attributes
            for att, val in atts.items():
                if isinstance(val, str):
                    atts[att] = val.format(**self.atts)

            # Set
            set_atts(data, atts)

        return data


    def unpack(self, pdata, rescale=True, format=1, firstdims=None, id=None,
            atts=None, format_atts=True):
        """Unpack data along space, reshape, and optionally unnormalize and remean.

        Input is sub_space:other, output is other:split_space.

        Params
        ------
        pdata: 2D array
            Packed data with space dim as first.
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
        # - create variable
        data = self.create_array(firstdims=firstdims, format=format,
            pshape=pdata.shape, id=id, atts=atts, format_atts=format_atts)
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

class Simple1DPacker(object):
    """Packer pure numpy fortran arrays along the first axis"""

    def __init__(self, data, missing_value, valid=None, **kwargs):

        if valid is None:
            self.valid = ~N.isclose(data, missing_value, **kwargs).any(axis=0)
        else:
            self.valid = valid
        self.compressible = not valid.all()
        self.packed_data = self.pack(data)
        self.missing_value = missing_value

    def pack(self, data, copy=False):
        if not self.compressible:
            if copy:
                data = data.copy()
            return data
        return N.asfortranarray(data[self.valid])

    def unpack(self, pdata, copy=False):
        if not self.compressible:
            if copy:
                pdata = pdata.copy()
            return pdata
        data = N.zeros(self.valid.shape[:1] + pdata.shape[1:])
        data += self.missing_value
        data[self.valid] = pdata
        return N.asfortranarray(data)

