# All tests in the right order
TESTS=$(shell python -c "from __init__ import print_test_order;print_test_order('prefix')")

# Other options
PYTHONOPTIONS=-W ignore

# Default target
all: 
	nosetests $(patsubst %,%.py,$(TESTS))

# Run a test
test_%: test_%.py
	nosetests $<
