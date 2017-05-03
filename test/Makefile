# All tests in the right order
TESTS=$(shell python -c "from __init__ import print_test_order;print_test_order('prefix')")

# Other options
PYTHONOPTIONS=-W ignore

# Default target
all: 
	export MPLBACKEND=agg
	export UVCDAT_ANONYMOUS_LOG=no
	nosetests -v $(patsubst %,%.py,$(TESTS))

# Run a test
test_%: test_%.py
	export MPLBACKEND=agg
	export UVCDAT_ANONYMOUS_LOG=no
	nosetests $<
