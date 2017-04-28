import os
from StringIO import StringIO
import contextlib
from sphinx.util.console import bold
from vcmq import checkdir, opt2rst
from sonat.cui import main as sonat_main

@contextlib.contextmanager
def capture():
    import sys
    from cStringIO import StringIO
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out=[StringIO(), StringIO()]
        sys.stdout,sys.stderr = out
        yield out
    finally:
        sys.stdout,sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()



def write_usages(app):
    for name in app.builder.status_iterator(
            app.config.sonathelps_commands,
            bold("generating helps... "),
            length=len(app.config.sonathelps_commands),
            ):

        # Get help string
        opts = app.config.sonathelps_commands[name]
        with capture() as out:
            try:
                sonat_main(opts)
            except:
                pass

        # Write help as rst
        help_file =app.config.sonathelps_help_pat.format(name)
        if help_file:
            help_file = os.path.join(app.builder.srcdir, help_file)
            checkdir(help_file, asfile=True)
            prog = ' '.join(['sonat'] + opts[:-1])
            f = open(help_file, 'w')
            f.write(opt2rst(out[0], prog=prog))
            f.close()

        # Write usage
        usage_file = app.config.sonathelps_usage_pat.format(name)
        if usage_file:
            usage_file = os.path.join(app.builder.srcdir, usage_file)
            checkdir(usage_file, asfile=True)
            f = open(usage_file, 'w')
            f.write(out[0])
            f.close()




def setup(app):
    app.add_config_value('sonathelps_commands', {}, '')
    app.add_config_value('sonathelps_help_pat', 'generated/help.{}.txt', '')
    app.add_config_value('sonathelps_usage_pat', 'generated/usage.{}.txt', '')

    app.connect('builder-inited', write_usages)

    return {'version': '0.1'}


