from __future__ import absolute_import

from .ppl_decl cimport *
from .linear_algebra cimport *

cdef _wrap_Generator(PPL_Generator generator)

cdef _wrap_Generator_System(PPL_Generator_System generator_system)

cdef PPL_GeneratorType_str(PPL_GeneratorType t)

cdef class Generator(object):
    cdef PPL_Generator *thisptr

cdef class Generator_System(object):
    cdef PPL_Generator_System *thisptr

cdef class Poly_Gen_Relation(object):
    cdef PPL_Poly_Gen_Relation *thisptr
