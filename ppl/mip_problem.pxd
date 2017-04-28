from .linear_algebra cimport *
from .generator cimport *
from .constraint cimport *

cdef class MIP_Problem(object):
    cdef PPL_MIP_Problem *thisptr
