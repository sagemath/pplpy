from ppl_decl cimport * 

cdef class _mutable_or_immutable(object):
    cdef bint _is_mutable

cdef class Variable(object):
    cdef PPL_Variable *thisptr

cdef class Linear_Expression(object):
    cdef PPL_Linear_Expression *thisptr

cdef class Poly_Gen_Relation(object):
    cdef PPL_Poly_Gen_Relation *thisptr

cdef class Poly_Con_Relation(object):
    cdef PPL_Poly_Con_Relation *thisptr
    
cdef _wrap_Constraint(PPL_Constraint constraint)

cdef _wrap_Constraint_System(PPL_Constraint_System constraint_system)
