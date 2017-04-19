import os
from StringIO import StringIO
from sphinx.util.console import bold
from vcmq import checkdir, cfg2rst
from sonat.config import SONAT_CFGM

def write_content(app):
    app.info(bold('writing default config to rst... '), nonl=True)
    config_rstfile = os.path.join(app.env.srcdir, app.config.sonatconfig_content_file)
    checkdir(config_rstfile)
    s = StringIO()
    SONAT_CFGM.defaults().write(s)
    s = s.getvalue()
    s = s.replace(os.environ['USER'], '${USER}')
    f = open(config_rstfile,  'w')
    f.write(s)
    f.close()
    app.info('done into '+os.path.basename(config_rstfile))

def write_options(app):
    app.info(bold('writing config options to rst... '), nonl=True)
    config_directives = os.path.join(app.env.srcdir, app.config.sonatconfig_options_file)
    checkdir(config_directives)
    f = open(config_directives, 'w')
    f.write(SONAT_CFGM.get_rst(mode='specs'))
    f.close()
    app.info('done into '+os.path.basename(config_directives))

def setup(app):
    app.add_config_value('sonatconfig_content_file', 'generated/config.defaults.txt', '')
    app.add_config_value('sonatconfig_options_file', 'generated/config.options.txt', '')

    app.connect('builder-inited', write_content)
    app.connect('builder-inited', write_options)

    return {'version': '0.1'}


