# Minimal makefile for Sphinx documentation
# NOTE: we removed the shortcut $(O) for passing options to
# sphinx-build as this is undocummented and sometimes
# confusing.

# You can set these variables from the command line.
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = pplpy
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)
