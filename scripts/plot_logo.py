import os
import pylab as P
import numpy as N
import cdms2
from vcmq import add_shadow, checkdir

dpi = 100
#color = 'C0'
color = '#114568'
xsize = 800
figpat = '../doc/source/images/sonat-{xsize}-{bgcolor}{shadow}.svg'

aspect = 1 / 4.5

xisize = 1. * xsize / dpi
yisize = xisize * aspect

sonat_size = 145. * xisize / 8.
o_size = 50 * xisize / 8.

figfiles = []
P.figure(figsize=(xisize, yisize), dpi=dpi)
kw = dict(ha='center', va='center', color=color, weight='bold')
ts = P.figtext(.5, .41, 'SONAT', size=sonat_size, **kw)
to = P.figtext(.302, .49, 'O', size=o_size, **kw)
shadow = ''

bgcolor = 'transp'
figfiles.append(figpat.format(**locals()))
checkdir(figfiles[-1], asfile=True)
P.savefig(figfiles[-1], transparent=True)

bgcolor = 'white'
figfiles.append(figpat.format(**locals()))
P.savefig(figfiles[-1])

bgcolor = 'blue'
figfiles.append(figpat.format(**locals()))
P.savefig(figfiles[-1], facecolor="#2882B9")

shadow = True
add_shadow(ts, width=9, alpha=.5, xoffset=6, yoffset=-6)
add_shadow(to, width=9, alpha=.5, xoffset=6, yoffset=-6)
P.axis('off')
shadow='-shadow'

bgcolor = 'transp'
figfiles.append(figpat.format(**locals()))
P.savefig(figfiles[-1], transparent=True)

bgcolor = 'white'
figfiles.append(figpat.format(**locals()))
P.savefig(figfiles[-1])

bgcolor = 'blue'
figfiles.append(figpat.format(**locals()))
P.savefig(figfiles[-1], facecolor="#2882B9")
#P.show()


for path in figfiles:
    for newxsize in 150, 200, 600, 800:
        newysize = newxsize * aspect
        newpath = path.replace(str(xsize), str(newxsize)).replace('.svg', '.png')
        option = '-transparent white' if 'transp' in path else ''
        cmd = 'convert {option} -resize {newxsize}x{newysize} {path} {newpath}'.format(**locals())
        os.system(cmd)
        print cmd

