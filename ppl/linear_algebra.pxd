from ppl_decl cimport *

cdef class Variable(object):
    cdef PPL_Variable *thisptr

cdef class Linear_Expression(object):
    cdef PPL_Linear_Expression *thisptr
