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

    def __init__(self, dataset, weights=None, norms=None,
        keep_invalids=False, minvalid=None, clean_weights=True,
        logger=None, loglevel=None, zerofill=False, **kwargs):

        # Logger
        Logger.__init__(self, logger=logger, loglevel=loglevel,
            **dict_filter(kwargs, 'log_'))

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
        self.nt = None
        weights = self.remap(weights, reshape=True)
        norms = self.remap(norms, reshape=True)
        if self.ndataset==1 and norms[0] is None: norms = [False]
        self._invalids = []
        self.masked = False

        # Loop on datasets
        for idata,data in enumerate(dataset):

            # Create the Data instance and pack array
            dd = Data(data, norm=norms[idata], weights=weights[idata],
                keep_invalids=keep_invalids, minvalid=minvalid, clean_weights=clean_weights,
                zerofill=zerofill)
            self.data.append(dd)
            self._invalids.append(dd.invalids)
            self.masked |= dd.masked

            # Check nt
            if self.nt is None:
                self.nt = dd.nt
            elif self.nt != dd.nt:
                self.error('Time dimension of variable %i must have length %i (not %i)'%(idata, self.nt, dd.nt))

        # Merge
        self.stacked_data = npy.asfortranarray(npy.vstack([d.packed_data for d in self.data]))
        self.splits = npy.cumsum([d.packed_data.shape[0] for d in self.data[:-1]])
        self.stacked_weights = npy.hstack([d.packed_weights for d in self.data])
        self.ns = self.stacked_data.shape[0]
        self.ntv = (self.stacked_data!=default_missing_value).any(axis=0).sum()

    def get_norms(self, idata=None):
        """Get :attr:`norms` for one or all input variables"""
        if idata is None:
            return [d.norm for d in self.dataset.data]
        return self[idata].norm
    norms = property(get_norms, doc="Normalization factors")

    def get_invalids(self, stacked=True):
        if self._invalids[0] is None: return
        if not stacked: return self._invalids
        return self.restack(self._invalids, scale=False)
    invalids = property(fget=get_invalids, doc='Final mask of stacked data')

    def restack(self, dataset, scale=True):
        """Stack new variables as a fortran array

        It has the opposite effect of :meth:`unstack`.

        :Params:

            - **dataset**: Argument in the same form as initialization data.
            - **scale**, optional: Scale the variable (mean and norm), and optionally
              remove spatial mean if == 2.

        :Seel also: :meth:`Data.repack`
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

        # Check time length (first axis)
        nt1 = dataset[0].size/self[0].nstot
        if len(dataset)>1:
            for i, d in enumerate(dataset[1:]):
                i += 1
                nt2 = d.size/self[i].nstot
                if packs[0].ndim!= packs[i].ndim:
                    if packs[0].ndim==1:
                        SpanlibError('Variable %i must not have a time dimension '% (i+1, ))
                    else:
                        SpanlibError('Variable %i must have a time dimension '% (i+1, ))
                elif nt1!=nt2:
                    self.error('Time length of variable %i (%i) different from that of first variable (%i)'
                        % (i+1, nt2, nt1))

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

    def fill_invalids(self, dataref, datafill, raw=False,  copy=False, unmap=True,
        missing=False):
        """Fill ``dataref`` with ``datafill`` at registered invalid points
        or current missing points


        :Params:

            - **dataref**: Reference dataset than must be filled.
            - **datafill**: Dataset used to fill reference.
            - **raw**, optional: Input data (and mask through missing) are aleardy packed.
            - **missing**, optional: If True, gaps are defined by ``dataref``
              missing values. If False, gaps are those of original data. If an array,
              it is used as the mask to define gaps.
        """
#        if self.invalids is None:
#            return dataref.clone() if copy else dataref # FIXME:

        # Re stack to have raw data
        if not raw:
            if dataref!='stacked_data':
                dataref = self.restack(dataref)
            datafill = self.restack(datafill)
            copy = False
        if dataref=='stacked_data':
             dataref = self.stacked_data

        # Check raw shape
        if dataref.shape!=datafill.shape:
            self.error("Can't replace with raw data because of wrong shape (%s!=%s)"%
                (data.shape, datafill.shape))

        # Copy for safety
        if copy: dataref = dataref.copy()

        # Put it at invalid points
        if missing is not False and missing is not True:
            if raw:
                mask = missing
            else:
                mask = self.restack(missing)
        elif missing:
            mask = npy.ma.masked_values(dataref, default_missing_value, shrink=False).mask
        else:
            mask = self.invalids
        dataref[:] = npy.where(mask, datafill, dataref)
        del mask

        # Unstack ?
        if int(raw)>0:
            data = npy.asfortranarray(dataref)
        else:
            data = self.unstack(dataref)
            if unmap: data = self.unmap(data)

        return data


    def replace_invalids(self, datafill, raw=False, copy=False):
        if self.invalids is None: return

        # Re stack to have raw data if needed
        if not raw:
            datafill = self.restack(datafill)

        # Fill
        return self.fill_invalids(self.stacked_data, datafill, raw=True, copy=False)




#    def get_time(self, idata=0):
#        """Get time axis"""
#        return self.data[idata].get_time()



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

    def get_time(self, nt=None, offset=0, idata=0):
        """Get the time axis of an input variable

        If CDAT is not used, length of axis is returned.
        """
        return self[idata].get_time(nt, offset)



