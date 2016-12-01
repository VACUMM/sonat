################################################################################
NAME=sonat
CURDIRNAME=$(shell basename $(CURDIR))
PYTHON_PACKAGE_NAME=$(NAME)
PYTHON_VERSION=$(shell  python -c 'import sys; print "%d.%d"%sys.version_info[:2]')
VERSION=$(shell python -c 'import setup; print setup.version')
RELEASE=${VERSION}
DATE=$(shell python -c 'import setup; print setup.date')
SDIST_FILE_NAME=$(NAME)-$(VERSION).tar.gz
SDIST_FILE_PATH=dist/$(SDIST_FILE_NAME)
TEST_INST_DIR=$(CURDIR)/test/install
TEST_SETUP_INST_DIR=$(TEST_INST_DIR)/setup-inst
#EXCLUDES=--exclude var --exclude *.pid --exclude *.log* --exclude *.out --exclude *.pyc
EXCLUDES=--exclude var

################################################################################
.PHONY: lib doc html pdf

################################################################################
all: help
help:
	@echo ""
	@echo "$(NAME) $(VERSION).$(RELEASE) ($(DATE))"
	@echo ""
	@echo "Usage: make [<target1> [<target2>...]]"
	@echo ""
	@echo "Available make targets:"
	@echo ""
	@echo "  Build distribution:"
	@echo "    sdist                     build python source distribution $(SDIST_FILE_PATH)"
	@echo "    dist                      call all build targets"
	@echo ""
	@echo "  Test:"
	@echo "    test-unittests            perform all unit tests using nosetests'
	@echo "    test-install              test python setup in $(TEST_SETUP_INST_DIR)"
	@echo "    test                      call all test targets"
	@echo ""
	@echo "  Clean"
	@echo "    clean-doc                 clean documentations"
	@echo "    clean-lib                 clean source package"
	@echo "    clean-build               clean builds"
	@echo "    clean-test                clean test installations"
	@echo "    clean-all                 call all clean and uninstall targets"
	@echo ""
	@echo "  Development"
	@echo "    lib                       local compilation of modules and extensions"
	@echo "    doc                       generate documentations"
	@echo "    html                      generate html documentation"
	@echo "    pdf                       generate pdf documentation"
	@echo "    safedoc                   generate documentations and all their dependencies"
	@echo "    install                   prepare for local use (build libs, fix permissions)"
	@echo "    uninstall                 clean local installations"
	@echo "    arch                      clean all and create source archive in parent directory"
	@echo ""

################################################################################
# DVELOPMENT
################################################################################

doc:
	cd doc && make

html:
	cd doc && make html

pdf:
	cd doc && make pdf

safedoc:
	make doc

lib:
	python setup.py build_ext --inplace --force

# install: lib doc
install: lib

uninstall: clean-lib clean-doc

arch: clean-all
	cd .. && tar cvjf $(CURDIRNAME).tbz $(CURDIRNAME) $(EXCLUDES)

################################################################################
# BUILD DIST
################################################################################

sdist: clean-partial install
	python setup.py sdist

dist: sdist

################################################################################
# TEST
################################################################################

test-install: clean-build clean-test-install
	python setup.py install -O1 --prefix=$(TEST_SETUP_INST_DIR)
	PYTHONPATH=$(TEST_SETUP_INST_DIR)/lib/python$(PYTHON_VERSION)/site-packages python -c "import "$(PYTHON_PACKAGE_NAME)"; "$(PYTHON_PACKAGE_NAME)".info()"

test-unittests:
	nosetests

test: test-unittests test-install

################################################################################
# CLEAN
################################################################################

clean-doc:
	-cd doc && make clean

clean-lib:
	-find lib/python -name '*.py[co]' -delete

clean-build:
	-rm -rf build MANIFEST setup.pyc

clean-dist:
	-rm -rf dist

clean-test-install:
	-rm -rf $(TEST_SETUP_INST_DIR)

clean-test: clean-test-install
	-rmdir $(TEST_INST_DIR)

clean-partial: clean-build clean-lib

clean: clean-partial clean-test

clean-all: uninstall clean clean-dist


