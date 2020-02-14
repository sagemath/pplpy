from __future__ import absolute_import

from .ppl_decl cimport *

cdef class Variable(object):
    cdef PPL_Variable *thisptr

cdef class Variables_Set(object):
    cdef PPL_Variables_Set *thisptr

cdef class Linear_Expression(object):
    cdef PPL_Linear_Expression *thisptr

cdef PPL_Coefficient PPL_Coefficient_from_pyobject(c) except *
