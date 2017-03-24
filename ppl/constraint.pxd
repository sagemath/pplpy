from .wrappers cimport *

cdef class Constraint(object):
    cdef PPL_Constraint *thisptr

cdef class Constraint_System(_mutable_or_immutable):
    cdef PPL_Constraint_System *thisptr

cdef class Constraint_System_iterator(object):
    cdef Constraint_System cs
    cdef cs_iterator_ptr csi_ptr
    