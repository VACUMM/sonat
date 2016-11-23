"""Test script for module :mod:`pyarm.misc`"""

from vcmq import adatetime, comptime

from util import (THISDIR, NCPAT_MANGA, NCFILE_MANGA0, NCFILE_MANGA1,
    NCGLOB_MANGA, NCPATGLOB_MANGA,
    assert_allclose, assert_raises, N, cdms2, cdtime)

from pyarm.misc import (scan_format_string, list_files_from_pattern,
    DatePat2GlobFormatter, ncfiles_time_indices)

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


if __name__=='__main__':
    test_misc_scan_format_string()
    test_misc_list_files_from_pattern()
    test_misc_datepat2globformatter()
    test_misc_ncfiles_time_indices()
