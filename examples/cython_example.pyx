r"""
Example of using pplpy from Cython.

Building this file needs the include paths of pplpy, PPL, GMP, cysignals and
gmpy2; see ``meson.build`` in this directory. ``meson test`` compiles it and
checks the results with ``test_cython_example.py``.
"""
#*****************************************************************************
#       Copyright (C) 2026 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from ppl.linear_algebra cimport Variable
from ppl.constraint cimport Constraint_System
from ppl.mip_problem cimport MIP_Problem
from ppl.polyhedron cimport C_Polyhedron


def mip_problem():
    r"""
    Return the problem of maximizing x + y with x, y >= 0 and 3x + 5y <= 10.
    """
    cdef Variable x = Variable(0)
    cdef Variable y = Variable(1)
    cdef Constraint_System cs = Constraint_System()
    cs.insert(x >= 0)
    cs.insert(y >= 0)
    cs.insert(3 * x + 5 * y <= 10)
    return MIP_Problem(2, cs, x + y)


def polyhedron():
    r"""
    Return the polyhedron from the example in the README.
    """
    cdef Variable x = Variable(0)
    cdef Variable y = Variable(1)
    cdef Variable z = Variable(2)
    cdef Constraint_System cs = Constraint_System()
    cs.insert(x >= 0)
    cs.insert(y >= 0)
    cs.insert(x + y + z == 1)
    return C_Polyhedron(cs)


def space_dimension(C_Polyhedron poly):
    return poly.thisptr.space_dimension()
