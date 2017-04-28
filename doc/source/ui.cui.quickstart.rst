Quick start
===========

.. highlight:: bash


Make sure to have a minimal command line help::

    $ sonat --short-help

Print SONAT info::

    $ sonat info

Try to open the web help with some a argument::

    $ sonat help NcObsPlatform

Create a directory where to work::

    $ mkdir work
    $ cd work

Copy the sample config file::

    $ cp $(sonat info test_dir)/sonat.cfg .

List sample data::

    $ ls $(sonat info test_dir)

Plot the available observations platforms setting the size of markers::

    $ sonat obs plot --obs-lots-size=18

Results are in the directory :file:`CUI/` by default.
Open the html file::

    $ xdg-open CUI/OBS/sonat.obs.html

Generate a pseudo-ensemble::

    $ sonat ens gen_pseudo

Plot ensemble diagnostics with more logging and with observation locations::

    $ sonat --logger-level=debug plot_diags --add-obs
    $ xdg-open CUI/ENS/sonat.ens.html

Make an ARM analysis::

    $ sonat arm analysis
    $ xdg-open CUI/ARM/sonat.arm.html
    
Make an ARM sensitivity analysis::

    $ sonat arm sa
    $ xdg-open CUI/ARM/sonat.armsa.html

