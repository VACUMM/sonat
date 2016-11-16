"""Stacking module

Inspired from the spanlib library.
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

from vcmq import dict_filter

class Dataset(_Base__):
    """Class to handle one variable or a list of variables

    This fonction concatenates several dataset that have the
    same time axis. It is useful for analysing for example
    several variables at the same time.
    It takes into account weights, masks and axes.

    :Params:

        - *dataset*: A list of data objects.
            They must all have the same time length.

    :Options:

        - *weights*: Associated weights.
    """

    def __init__(self, dataset, norms=None, **kwargs):

        # Logger
        _Base_.__init__(self, logger=logger, **dict_filter(kwargs, 'logger_'))

        # Input shape
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            self.map = len(dataset)
        else:
            dataset = [dataset]
            self.map = 0
        self.ndataset = self.nd = len(dataset)
        self.dataset = dataset

        # Other inits
        self.data = []
        self.nr = None
        norms = self.remap(norms, reshape=True)
        if self.ndataset==1 and norms[0] is None: norms = [False]
        self.masked = False

        # Loop on datasets
        for idata,data in enumerate(dataset):

            # Create the Data instance and pack array
            dd = Data(data, norm=norms[idata],
                keep_invalids=keep_invalids, minvalid=minvalid,
                zerofill=zerofill)
            self.data.append(dd)
            self.masked |= dd.masked

            # Check nr
            if self.nr is None:
                self.nr = dd.nr
            elif self.nr != dd.nr:
                self.error('Record dimension of variable %i must have length %i (not %i)'%(idata, self.nr, dd.nr))

        # Merge
        self.stacked_data = npy.asfortranarray(npy.vstack([d.packed_data for d in self.data]))
        self.splits = npy.cumsum([d.packed_data.shape[0] for d in self.data[:-1]])
        self.ns = self.stacked_data.shape[0]
        self.nrv = (self.stacked_data!=default_missing_value).any(axis=0).sum()

    def get_norms(self, idata=None):
        """Get :attr:`norms` for one or all input variables"""
        if idata is None:
            return [d.norm for d in self.dataset.data]
        return self[idata].norm
    norms = property(get_norms, doc="Normalization factors")

    def restack(self, dataset, scale=True):
        """Stack new variables as a fortran array

        It has the opposite effect of :meth:`unstack`.

        :Params:

            - **dataset**: Argument in the same form as initialization data.
            - **scale**, optional: Scale the variable (mean and norm), and optionally
              remove spatial mean if == 2.

        :Seel also: :meth:`Pack.repack`
        """

        # Check input data
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            dmap = len(dataset)
        else:
            dataset = [dataset]
            dmap = 0
        if  len(dataset)!=self.ndataset:
            self.error('You must provide %i variable(s) to stack'%len(dataset))

        # Pack
        packs = [self[idata].repack(data, scale=scale, force2d=True)
            for idata, data in enumerate(dataset)]

        # Check record length (first axis)
        nr1 = dataset[0].size/self[0].nstot
        if len(dataset)>1:
            for i, d in enumerate(dataset[1:]):
                i += 1
                nr2 = d.size/self[i].nstot
                if packs[0].ndim!= packs[i].ndim:
                    if packs[0].ndim==1:
                        SpanlibError('Variable %i must not have a time dimension '% (i+1, ))
                    else:
                        SpanlibError('Variable %i must have a time dimension '% (i+1, ))
                elif nr1!=nr2:
                    self.error('Time length of variable %i (%i) different from that of first variable (%i)'
                        % (i+1, nr2, nr1))

        # Stack
        stacker = npy.vstack if packs[0].ndim==2 else npy.hstack
        sdata = npy.asfortranarray(stacker(packs))

        return sdata


    def unstack(self, sdata, rescale=True, format=1, firstdims=None, firstaxes=None):
        """Unstack and unpack data

        It has the opposite effect of :meth:`restack`.


        :Params:

            - **sdata**: Stacked data (as computed by :meth:`__init__` or :meth:`restack`).
            - **rescale**, optional: Rescale the variable (mean and norm) and
              add spatial mean if == 2.
            - **format**, optional: Format the variable (see :meth:`Data.create_array`).
            - **firstaxes**, optional: First axis (see :meth:`Data.create_array`).


        :Seel also: :meth:`Data.unpack`
        """

        # Unstack
        spliter = npy.hsplit if sdata.ndim==1 else npy.vsplit
        packs = spliter(sdata, self.splits)

        # Unpack
        return [self[i].unpack(pdata, rescale=rescale, format=format,
            firstdims=firstdims, firstaxes=firstaxes)
            for i, pdata in enumerate(packs)]


    def __len__(self):
        return self.ndataset

    def __getitem__(self, i):
        return self.data[i]

    def unmap(self, out):
        """Remap out so as to be in the same form as input data

        It has the opposite effect of :meth:`remap`.
        """
        if self.map==0:
            return out[0]
        return out

    def remap(self, values, reshape=True, fill_value=None):
        """Makes sure that values is a list (or tuple) of length :attr:`ndataset`

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
                self.error('Wrong number of input items (%i instead of %i)'
                    %(len(values), n))
            values = broadcast(values, n, fill_value)
        return values

    def has_cdat(self, idata=None):
        """Check if there were at least one input array of CDAT type (:mod:`MV2`)

        :Sea also: :meth:`Data.has_cdat`
        """
        # Not at all
        if not has_cdat_support: return False

        # Single var
        if idata is not None:
            return self[idata].has_cdat()

        # At least one
        for d in self.data:
            if d.has_cdat(): return True
        return False

    def get_time(self, nr=None, offset=0, idata=0):
        """Get the time axis of an input variable

        If CDAT is not used, length of axis is returned.
        """
        return self[idata].get_time(nr, offset)



