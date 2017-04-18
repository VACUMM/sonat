import os
from StringIO import StringIO
from sphinx.util.console import bold
from vcmq import checkdir



def write_ncdumps(app):
    for name in app.builder.status_iterator(
            app.config.sonatncdumph_ncfiles,
            bold("generating list of ncdump -h... "),
            length=len(app.config.sonatncdumph_ncfiles),
            ):
        txtfile = app.config.sonatncdumph_filepat.format(name)
        txtfile = os.path.join(app.builder.srcdir, txtfile)
        checkdir(txtfile, asfile=True)
        ncfile = app.config.sonatncdumph_ncfiles[name]
        if not os.path.isabs(ncfile):
            ncfile = os.path.join(app.builder.srcdir, ncfile)
        os.system('ncdump -h {ncfile} > {txtfile}'.format(**locals()))



def setup(app):
    app.add_config_value('sonatncdumph_ncfiles', {}, '')
    app.add_config_value('sonatncdumph_filepat', 'generated/ncdump.{}.txt', '')

    app.connect('builder-inited', write_ncdumps)

    return {'version': '0.1'}


