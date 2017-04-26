from .linear_algebra cimport *
from .generator cimport *
from .constraint cimport *

cdef class MIP_Problem(object):
    cdef PPL_MIP_Problem *thisptr

cdef class MIP_Problem_constraints_iterator(object):
    cdef MIP_Problem pb
    cdef PPL_mip_iterator *mip_csi_ptr
