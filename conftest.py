r"""
Collects the doctests in pplpy.

The usual --doctest-modules flag of pytest doesn't work here since it doesn't
consider tests in .pyx files.
"""
#*****************************************************************************
#       Copyright (C) 2026 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import importlib

from _pytest.doctest import DoctestModule


class NamedDoctestModule(DoctestModule):
    r"""
    Collects doctests for a module that is imported by name and not by path.

    This allows us to load compiled Cython modules and also the installed
    __init__.py of ppl.

    Currently, this only considers modules that live directly in ppl/. If you
    add more structure to the pplpy package, the logic here needs to be
    expanded.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        module = ["ppl"]

        if self.path.name != "__init__.py":
            module.append(self.path.stem)

        self._module_name = ".".join(module)

    def _getobj(self):
        r"""
        Return the loaded module so pytest can extract doctests from its docstrings.
        """
        return importlib.import_module(self._module_name)


def pytest_collect_file(file_path, parent):
    r"""
    Return a doctest collector for the modules of the ``ppl`` package.

    pytest calls this hook for every file it visits. Files directly inside
    ``ppl/`` that are ``.pyx`` sources or ``__init__.py`` are mapped to their
    module name; for everything else, return ``None`` and leave the file to
    pytest's other collectors.
    """
    if file_path.relative_to(parent.config.rootpath).parts[0] == "ppl" and file_path.suffix in [".pyx", ".py"]:
        return NamedDoctestModule.from_parent(parent, path=file_path)
