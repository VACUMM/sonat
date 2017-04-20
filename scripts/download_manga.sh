#!/bin/bash

# MANGA
#    Start: 2006-01-01T00:00:00Z
#    End: 2015-01-15T00:00:00Z
#    DT: 1 hour

URL="http://tds1.ifremer.fr/thredds/dodsC/PREVIMER-MANGAE4000-MARS3DF1-FOR_FULL_TIME_SERIE"
DATE00=2014-01-01
DATE01=2014-01-15
DATE10=2014-01-16
DATE11=2014-02-01
I0=226
I1=292
J0=174
J1=204
IJSAMP=2
TSAMP=24
KSAMP=2
NCFILE0="../data/manga-$DATE00.nc"
NCFILE1="../data/manga-$DATE10.nc"

ISEL=$I0,$I1
JSEL=$J0,$J1
XYSEL="-d ni,$ISEL -d nj,$JSEL -d ni_u,$ISEL -d nj_u,$JSEL -d ni_v,$ISEL -d nj_v,$JSEL -d ni_f,$ISEL -d nj_f,$JSEL"
SAMP="-d time,,,$TSAMP -d ni,,,$IJSAMP -d nj,,,$IJSAMP -d ni_u,,,$IJSAMP -d nj_u,,,$IJSAMP -d ni_v,,,$IJSAMP -d nj_v,,,$IJSAMP -d ni_f,,,$IJSAMP -d nj_f,,,$IJSAMP -d level,,,$KSAMP"

set -x

ncks -O -d time,$DATE00,$DATE01 $XYSEL $URL tmp.nc || exit
ncks -O -4 $SAMP tmp.nc $NCFILE0

ncks -O -d time,$DATE10,$DATE11 $XYSEL $URL tmp.nc || exit
ncks -O -4 $SAMP tmp.nc $NCFILE1

rm -rf tmp.nc



