from .ppl_decl cimport *

cdef _wrap_Congruence(PPL_Congruence)

cdef class Congruence(object):
    cdef PPL_Congruence *thisptr

cdef class Congruence_System(object):
    cdef PPL_Congruence_System *thisptr
