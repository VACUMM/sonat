Introduction
============

SONAT is an open-source python framework intended to help designing a multivariate, multiplatform
ocean observation network using ocean model ensembles.

With SONAT, a given network is evaluated using the so-called ArM (Array Modes)
method :cite:`bene85,bene90,leheal09,lamoal16`.
The ArM method is a light-weight approach to network design,
based on the ability of a configuration to capture most
of the model error covariances.
These errors as simulated using ensemble methods.
An "good" network is supposed to maximise the impact of observations
in an assimilation system.
However, unlike the so-called OSSE (Observing System Simulation Experiments) approach,
no assmilation is needed to assess a network. 

In the current implementation of SONAT, the ensemble is not representative
of model errors, and instead of model variance.
We speek of a 'pseudo-ensemble' that has no time dimension.
The ArM method becomes a measure of the efficiency of the network
to sample the variability of the system, without any perspective
of data assimilation.
SONAT facilitates the generation of pseudo-ensemble from model outputs.

An observation network may be constituted of very
different platforms such as profiles, HF radars, satellite products
or glider data, bottom measurements,
SONAT has been designed to handle such platform
types in a simple and generic way.

Depending on the application, very different variable types
may have to considered at the same time for network designing.
SONAT offer a seamless integration of variables of different types.

SONAT comes with a library, a commandline user interface and a graphical
user interface, and is configurable and customisable.
