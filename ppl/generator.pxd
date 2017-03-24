from ppl_decl cimport *
from .wrappers cimport *

cdef _wrap_Generator(PPL_Generator generator)

cdef _wrap_Generator_System(PPL_Generator_System generator_system)

cdef PPL_GeneratorType_str(PPL_GeneratorType t)

cdef class Generator(object):
    cdef PPL_Generator *thisptr

cdef class Generator_System(_mutable_or_immutable):
    cdef PPL_Generator_System *thisptr

cdef class Generator_System_iterator(object):
    cdef Generator_System gs
    cdef gs_iterator_ptr gsi_ptr