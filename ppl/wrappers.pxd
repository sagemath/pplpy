from ppl_decl cimport *

cdef PPL_GeneratorType_str(PPL_GeneratorType t)

cdef _wrap_Generator(PPL_Generator generator)

cdef _wrap_Constraint(PPL_Constraint constraint)

cdef class _mutable_or_immutable(object):
    cdef bint _is_mutable

cdef class Variable(object):
    cdef PPL_Variable *thisptr

cdef class Linear_Expression(object):
    cdef PPL_Linear_Expression *thisptr

cdef class Generator(object):
    cdef PPL_Generator *thisptr

cdef class Generator_System(_mutable_or_immutable):
    cdef PPL_Generator_System *thisptr

cdef class Generator_System_iterator(object):
    cdef Generator_System gs
    cdef gs_iterator_ptr gsi_ptr

cdef class Constraint(object):
    cdef PPL_Constraint *thisptr

cdef class Constraint_System(_mutable_or_immutable):
    cdef PPL_Constraint_System *thisptr

cdef class Constraint_System_iterator(object):
    cdef Constraint_System cs
    cdef cs_iterator_ptr csi_ptr

cdef class Polyhedron(_mutable_or_immutable):
   cdef PPL_Polyhedron *thisptr
   cdef _relation_with_generator(Polyhedron self, Generator g)
   cdef _relation_with_constraint(Polyhedron self, Constraint c)

cdef class C_Polyhedron(Polyhedron):
    pass

cdef class NNC_Polyhedron(Polyhedron):
    pass

cdef class Poly_Gen_Relation(object):
    cdef PPL_Poly_Gen_Relation *thisptr

cdef class Poly_Con_Relation(object):
    cdef PPL_Poly_Con_Relation *thisptr
