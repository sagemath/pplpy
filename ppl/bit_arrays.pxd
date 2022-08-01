from .ppl_decl cimport *

cdef class Bit_Row(object):
    cdef PPL_Bit_Row *thisptr

cdef class Bit_Matrix(object):
    cdef PPL_Bit_Matrix *thisptr
