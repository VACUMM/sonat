import os
if not os.environ.get('DISPLAY', ''): # detect batch mode
    from matplotlib import use
    use('Agg', warn=False)
import sys
import MV2, cdms2
from vcmq import (lindates, create_lon, create_lat, create_dep,
    create_time, N, rotate_grid, set_grid, create_axis)

THIS_DIR = os.path.dirname(__file__)

# Try local import first
LIB_DIR = os.path.abspath(os.path.join(THIS_DIR, '..', 'sonat'))
if os.path.exists(os.path.join(LIB_DIR, '__init__.py')):
    sys.path.insert(0, os.path.dirname(LIB_DIR))
    try:
        import sonat._fcore
    except:
        del sys.path[0]
        if os.environ.get('READTHEDOCS', ''):
            del sys.path[:2]
from sonat import get_logger, get_data_dir, LOGGER


assert_allclose = N.testing.assert_allclose
assert_raises = N.testing.assert_raises

DATA_DIR = get_data_dir()

NC_ENS = os.path.abspath(os.path.join(DATA_DIR, 'ens.nc'))
NCPAT_MANGA = os.path.abspath(os.path.join(DATA_DIR, 'manga-{date:%Y-%m-%d}.nc'))
NCFILE_MANGA0 = os.path.abspath(os.path.join(DATA_DIR, 'manga-2014-01-01.nc'))
NCFILE_MANGA1 = os.path.abspath(os.path.join(DATA_DIR, 'manga-2014-01-16.nc'))
NCGLOB_MANGA = os.path.abspath(os.path.join(DATA_DIR, 'manga-*-[01][0-9]-??.nc'))
NCPATGLOB_MANGA = os.path.abspath(os.path.join(DATA_DIR, 'manga-*-{date:%m-%d}.nc'))
NCFILE_OBS_SURF = os.path.abspath(os.path.join(DATA_DIR, 'obs.surf.nc'))
NCFILE_OBS_HFRADARS = os.path.abspath(os.path.join(DATA_DIR, 'obs.hfradars.nc'))
NCFILE_OBS_PROFILES = os.path.abspath(os.path.join(DATA_DIR, 'obs.profiles.nc'))
NCFILE_OBS_SATSST = os.path.abspath(os.path.join(DATA_DIR, 'obs.satsst.nc'))
NCFILE_BATHY = os.path.abspath(os.path.join(DATA_DIR, 'bathy.brittany.nc'))
NCFILE_BATHY_VARNAME = 'elevation'
NCFILE_BATHY_POSITIVE_UP = True
CFGFILE = os.path.abspath(os.path.join(THIS_DIR, 'sonat.cfg'))


if LOGGER is None:
    LOGGER_LEVEL_TESTS = os.environ.get("SONAT_LOGGER_LEVEL_TESTS", "error")
    LOGGER = get_logger(level=LOGGER_LEVEL_TESTS)

CACHE = {}

def get_bathy():
    """Positive up"""
    if 'bathy' not in CACHE:
        f = cdms2.open(NCFILE_BATHY)
        CACHE['bathy'] = f(NCFILE_BATHY_VARNAME)
        f.close()
        if not NCFILE_BATHY_POSITIVE_UP:
            CACHE['bathy'][:] *= -1
    return CACHE['bathy']

def create_mv2_gridder_xyzt(nx=8, ny=7, nz=6, nt=5, xmin=-6., xmax=-3, ymin=46,
        ymax=48, zmin=-200, zmax=0, tmin='2016', tmax='2016-02',
        tunits='days since 2016-01-01',
        rotate=0):
    """Create a MV2 array on a grid

    Return
    ------
    MV2.array
    """

    # Axes
    shape = ()
    axes = []
    if nt!=0:
        time = create_time(lindates(tmin, tmax, nt), tunits)
        axes.append(time)
        shape += nt,
    if nz!=0:
        dep = create_dep(N.linspace(zmin, zmax, nz))
        axes.append(dep)
        shape += nz,
    if ny!=0:
        lat = create_lat(N.linspace(ymin, ymax, ny))
        axes.append(lat)
        shape += ny,
    if nx!=0:
        lon = create_lon(N.linspace(xmin, xmax, nx))
        axes.append(lon)
        shape += nx,

    # Array
    data = MV2.array(N.arange(N.multiply.reduce(shape)).reshape(shape), copy=False,
        axes=axes, id='temp', dtype='d')

    # Rotate grid
    if rotate:
        grid = data.getGrid()
        if grid is not None:
            grid = rotate_grid(grid, rotate)
            set_grid(data, grid)

    return data

def create_mv2_scattered_xyzt(np=10, nz=6, nt=5, xmin=-6., xmax=-3, ymin=46,
        ymax=48, zmin=-200, zmax=0, tmin='2016', tmax='2016-02',
        tunits='days since 2016-01-01'):
    """Create a VM2 array of scattered data

    Return
    ------
    array: longitudes
    array: latitude
    MV2.array: data
    """

     # Axes
    shape = ()
    axes = []
    if nt!=0:
        time = create_time(lindates(tmin, tmax, nt), tunits)
        shape += nt,
        axes.append(time)
    if nz!=0:
        dep = create_dep(N.linspace(zmin, zmax, nz))
        axes.append(dep)
        shape += nz,
    shape += np,
    axes.append(create_axis((np, )))

    # Array
    data = MV2.array(N.arange(N.multiply.reduce(shape)).reshape(shape), copy=False,
        axes=axes, id='temp', dtype='d')

    # Positiions
    lons = N.linspace(xmin, xmax, np)
    lats = N.linspace(ymin, ymax, np)

    return lons, lats, data
