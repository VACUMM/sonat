
# Force the setup.cfg
cat > setup.cfg << EOF
[config_fc]
f90flags=-O3 -fdefault-real-8 -fbacktrace

[build_ext]
libraries=sangoma_tools lapack blas
EOF

# Prevent X and CDAT asking for usage logging
export MPLBACKEND=agg
export UVCDAT_ANONYMOUS_LOG=no

# Standard installation
$PYTHON -m pip install --no-deps --ignore-installed .
