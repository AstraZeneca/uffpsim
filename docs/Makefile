# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= $(shell if [ -x ../venv/bin/sphinx-build ]; then echo ../venv/bin/sphinx-build; elif command -v sphinx-build >/dev/null 2>&1; then echo sphinx-build; else echo sphinx-build; fi)
SOURCEDIR     = .
BUILDDIR      = _build
PYTHONPATH    = $(abspath ..)

export PYTHONPATH

.PHONY: help html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man gettext changes linkcheck doctest coverage clean

clean:
	rm -rf $(BUILDDIR)/*

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

dirhtml:
	@$(SPHINXBUILD) -M dirhtml "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

singlehtml:
	@$(SPHINXBUILD) -M singlehtml "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

pickle:
	@$(SPHINXBUILD) -M pickle "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

json:
	@$(SPHINXBUILD) -M json "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

htmlhelp:
	@$(SPHINXBUILD) -M htmlhelp "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

qthelp:
	@$(SPHINXBUILD) -M qthelp "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

devhelp:
	@$(SPHINXBUILD) -M devhelp "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

epub:
	@$(SPHINXBUILD) -M epub "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

latex:
	@$(SPHINXBUILD) -M latex "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

latexpdf:
	@$(SPHINXBUILD) -M latexpdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

text:
	@$(SPHINXBUILD) -M text "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

man:
	@$(SPHINXBUILD) -M man "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

gettext:
	@$(SPHINXBUILD) -M gettext "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

changes:
	@$(SPHINXBUILD) -M changes "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

linkcheck:
	@$(SPHINXBUILD) -M linkcheck "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doctest:
	@$(SPHINXBUILD) -M doctest "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

coverage:
	@$(SPHINXBUILD) -M coverage "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
