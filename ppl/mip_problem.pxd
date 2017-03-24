from .wrappers cimport *
from .generator cimport *
from .constraint cimport * 

cdef class MIP_Problem(object):
    cdef PPL_MIP_Problem *thisptr

cdef class MIP_Problem_constraints_iterator(_mutable_or_immutable):
    cdef MIP_Problem pb
    cdef mip_cs_iterator_ptr mip_csi_ptr