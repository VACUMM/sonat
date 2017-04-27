#!/bin/bash
#
# Installs conda (https://conda.io/miniconda.html) and its SONAT package

# Where to install
CONDA_DIR=$HOME/miniconda2
# 
# Install conda
KERNEL=$(uname -s)
HARDWARE=$(uname -i)
CONDA_INSTALLER=Miniconda2-latest-$KERNEL-$HARDWARE.sh
test ! -f $CONDA_INSTALLER && wget https://repo.continuum.io/miniconda/$CONDA_INSTALLER
chmod +x $CONDA_INSTALLER
$CONDA_INSTALLER -p $CONDA_DIR -b || exit 1
export PATH=$CONDA_DIR/bin:$PATH

# Install SONAT
conda install -y -c vacumm/label/test -c uvcdat -c conda-forge sonat || exit 1

# Info to test it
echo "IN AN EMPTY DIRECTORY, TEST WITH:"
echo "$CONDA_DIR/bin/sonat test"