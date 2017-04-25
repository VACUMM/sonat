"""Test script for module :mod:`sonat.misc`"""

import os
from vcmq import adatetime, comptime

from util import (THIS_DIR, NCPAT_MANGA, NCFILE_MANGA0, NCFILE_MANGA1,
    NCGLOB_MANGA, NCPATGLOB_MANGA, NC_ENS,
    assert_allclose, assert_raises, N, cdms2, cdtime)

from sonat.misc import (scan_format_string, list_files_from_pattern,
    DatePat2GlobFormatter, ncfiles_time_indices, slice_gridded_var,
    dicttree_relpath)

def test_misc_scan_format_string():

    fields, props = scan_format_string('myfile-{date:%Y-%m-%d}-{0:d}.nc')
    assert 'date' in fields
    assert 0 in fields
    assert 'date' in props['with_time']
    assert 'date' in props['keyword']
    assert 0 in props['positional']


def test_misc_list_files_from_pattern():

    files = list_files_from_pattern(NCFILE_MANGA0)
    assert files == [NCFILE_MANGA0]

    files = list_files_from_pattern([NCFILE_MANGA0, NCFILE_MANGA1])
    assert files == [NCFILE_MANGA0, NCFILE_MANGA1]

    files = list_files_from_pattern(NCPAT_MANGA)
    assert files == [NCFILE_MANGA0, NCFILE_MANGA1]

    files = list_files_from_pattern(NCPAT_MANGA, time=('2014-01-05', '2014-01-17'))
    assert files == [NCFILE_MANGA1]

    files = list_files_from_pattern(NCPAT_MANGA, time=('2014-01-05', '2014-01-17'),
        dtfile=(14, 'days'))
    assert files == [NCFILE_MANGA0, NCFILE_MANGA1]

    files = list_files_from_pattern(NCPAT_MANGA, time=('2014-01-05', '2014-01-12'),
        dtfile=(14, 'days'))
    assert files == [NCFILE_MANGA0]

    files = list_files_from_pattern(NCGLOB_MANGA)
    assert files == [NCFILE_MANGA0, NCFILE_MANGA1]

    files = list_files_from_pattern(NCPATGLOB_MANGA, time=('2014-01-05', '2014-01-17'),
        dtfile=(14, 'days'))
    assert files == [NCFILE_MANGA0, NCFILE_MANGA1]

def test_misc_datepat2globformatter():

    f = DatePat2GlobFormatter()
    assert f.format('aaa-{date:%Y-%m}.nc', date=adatetime('2000-10'))== "aaa-2000-10.nc"
    assert f.format('aaa-{date:%Y-%m}.nc') == "aaa-[0-2][0-9][0-9][0-9]-[0-1][0-9].nc"

def test_misc_ncfiles_time_indices():

    ncfiles = [NCFILE_MANGA0, NCFILE_MANGA1]
    dates = [comptime('2013'), adatetime('2014-01-6 15:24'),
        '2014-01-6 11:24', '2014-01-18', '2014-01-18 01','2016']
    ncfdict, info = ncfiles_time_indices(ncfiles, dates, getinfo=True)

    assert ncfdict[NCFILE_MANGA0] == [5, 6]
    assert ncfdict[NCFILE_MANGA1] == [2]
    assert info['duplicates'] == [comptime('2014-01-18 01')]
    assert info['missed'] == [comptime('2013'), comptime('2016')]


def test_slice_gridded_var():
    # Get vars
    f = cdms2.open(NCFILE_MANGA0)
    temp = f('TEMP', time=slice(0, 3))
    f.close()
    f = cdms2.open(NC_ENS)
    tempe = f('temp', member=slice(0, 3))
    f.close()

    # Slice it
    t = slice_gridded_var(temp, lon=-4.)
    assert t.getOrder() == 'tzy'
    t = slice_gridded_var(temp, lat=47.5)
    assert t.getOrder() == 'tzx'
    t = slice_gridded_var(temp, time='2014-01-01 15.')
    assert t.getOrder() == 'zyx'
    t = slice_gridded_var(temp, time='2014-01-01 15.', lon=-4.)
    assert t.getOrder() == 'zy'
    t = slice_gridded_var(temp, lat=47.5, lon=-4.)
    assert t.getOrder() == 'tz'
    t = slice_gridded_var(temp, lon=[-4., -4.1])
    assert t.getOrder() == 'tzyx' and t.shape[-1] == 2
    s = slice_gridded_var(tempe, depth=-12)
    assert s.getOrder() == '-yx'
    assert not s.mask.all()


def test_misc_dicttree_relpath():

    assert (dicttree_relpath(os.path.abspath('../data/myfile.nc'), '.')==
        '../data/myfile.nc')

    print dicttree_relpath(os.path.abspath('../data/myfile.nc'), '.')

    dd = {'sec1': {
        'key1':'ENS/sub/file.f90',
        'dict':{
            'key':os.path.abspath('../data/myfile.nc'),
            'list':[os.path.abspath('../data/ens/fig1.nc'), 'fig2.nc'],
            }
        }
    }
    print dicttree_relpath(dd, 'ENS')

if __name__=='__main__':
#    test_misc_scan_format_string()
#    test_misc_list_files_from_pattern()
#    test_misc_datepat2globformatter()
#    test_misc_ncfiles_time_indices()
    test_slice_gridded_var()
#    test_misc_dicttree_relpath()
