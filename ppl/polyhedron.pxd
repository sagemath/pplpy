from .wrappers cimport *
from .generator cimport * 
from .constraint cimport * 

cdef class Polyhedron(_mutable_or_immutable):
   cdef PPL_Polyhedron *thisptr
   cdef _relation_with_generator(Polyhedron self, Generator g)
   cdef _relation_with_constraint(Polyhedron self, Constraint c)

cdef class C_Polyhedron(Polyhedron):
    pass

cdef class NNC_Polyhedron(Polyhedron):
    pass