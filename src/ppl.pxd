from libcpp cimport bool as cppbool

from mpz cimport *
from gmpxx cimport *

####################################################
# Cython does not support ctypedef within cppclass; Hack around this restriction:
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::Generator":
    ctypedef enum PPL_GeneratorType "Parma_Polyhedra_Library::Generator::Type":
        LINE, RAY, POINT, CLOSURE_POINT

cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::Constraint":
    ctypedef enum PPL_ConstraintType "Parma_Polyhedra_Library::Constraint::Type":
        EQUALITY, NONSTRICT_INEQUALITY, STRICT_INEQUALITY

cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::MIP_Problem":
    ctypedef enum PPL_MIP_Problem_Control_Parameter_Name:
        PRICING
    ctypedef enum PPL_MIP_Problem_Control_Parameter_Value:
        PRICING_STEEPEST_EDGE_FLOAT, PRICING_STEEPEST_EDGE_EXACT, PRICING_TEXTBOOK


####################################################
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":

    ctypedef size_t PPL_dimension_type  "Parma_Polyhedra_Library::dimension_type"
    ctypedef mpz_class PPL_Coefficient  "Parma_Polyhedra_Library::Coefficient"
    cdef cppclass PPL_Variable          "Parma_Polyhedra_Library::Variable"
    cdef cppclass PPL_Linear_Expression "Parma_Polyhedra_Library::Linear_Expression"
    cdef cppclass PPL_Generator         "Parma_Polyhedra_Library::Generator"
    cdef cppclass PPL_Generator_System  "Parma_Polyhedra_Library::Generator_System"
    cdef cppclass PPL_Constraint        "Parma_Polyhedra_Library::Constraint"
    cdef cppclass PPL_Constraint_System "Parma_Polyhedra_Library::Constraint_System"
    cdef cppclass PPL_Polyhedron        "Parma_Polyhedra_Library::Polyhedron"
    cdef cppclass PPL_C_Polyhedron      "Parma_Polyhedra_Library::C_Polyhedron"   (PPL_Polyhedron)
    cdef cppclass PPL_NNC_Polyhedron    "Parma_Polyhedra_Library::NNC_Polyhedron" (PPL_Polyhedron)
    cdef cppclass PPL_Poly_Gen_Relation "Parma_Polyhedra_Library::Poly_Gen_Relation"
    cdef cppclass PPL_Poly_Con_Relation "Parma_Polyhedra_Library::Poly_Con_Relation"
    cdef cppclass PPL_MIP_Problem       "Parma_Polyhedra_Library::MIP_Problem"

    cdef cppclass PPL_Variable:
        PPL_Variable(PPL_dimension_type i)
        PPL_dimension_type id()
        bint OK()
        PPL_dimension_type space_dimension()

    cdef cppclass PPL_Linear_Expression:
        PPL_Linear_Expression()
        PPL_Linear_Expression(PPL_Linear_Expression &e)
        PPL_Linear_Expression(PPL_Coefficient n)
        PPL_Linear_Expression(PPL_Variable v)
        PPL_dimension_type space_dimension()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient inhomogeneous_term()
        bint is_zero()
        bint all_homogeneous_terms_are_zero()
        void ascii_dump()
        bint OK()
        PPL_Linear_Expression operator+(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator-(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator*(PPL_Coefficient n)
        PPL_Constraint operator> (PPL_Linear_Expression& e)
        PPL_Constraint operator>=(PPL_Linear_Expression& e)
        PPL_Constraint operator==(PPL_Linear_Expression& e)
        PPL_Constraint operator<=(PPL_Linear_Expression& e)
        PPL_Constraint operator< (PPL_Linear_Expression& e)

    cdef cppclass PPL_Generator:
        PPL_Generator(PPL_Generator &g)
        # Cython does not support static members
        #PPL_Generator line(PPL_Linear_Expression &e)
        #PPL_Generator ray(PPL_Linear_Expression &e)
        #PPL_Generator point(PPL_Linear_Expression &e, PPL_Coefficient d)
        #PPL_Generator closure_point(PPL_Linear_Expression &e)
        PPL_dimension_type space_dimension()
        PPL_GeneratorType type()
        bint is_line()
        bint is_ray()
        bint is_line_or_ray()
        bint is_point()
        bint is_closure_point()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient divisor() except +
        bint is_equivalent_to(PPL_Generator &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Constraint:
        PPL_Constraint(PPL_Constraint &g)
        PPL_dimension_type space_dimension()
        PPL_ConstraintType type()
        bint is_equality()
        bint is_inequality()
        bint is_nonstrict_inequality()
        bint is_strict_inequality()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient inhomogeneous_term()
        bint is_tautological()
        bint is_inconsistent()
        bint is_equivalent_to(PPL_Constraint &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Generator_System:
        # This seems to not work in cython
        #cppclass PPL_const_iterator "const_iterator":
        #    PPL_Generator operator*()
        #    PPL_const_iterator operator++()
        #    bint operator==(PPL_const_iterator&)
        #    bint operator!=(PPL_const_iterator&)
        #PPL_const_iterator begin()
        #PPL_const_iterator end()
        PPL_Generator_System()
        PPL_Generator_System(PPL_Generator &g)
        PPL_Generator_System(PPL_Generator_System &gs)
        PPL_dimension_type space_dimension()
        void clear()
        void insert(PPL_Generator &g)
        bint empty()
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Constraint_System:
        # This seems to not work in cython
        #cppclass PPL_const_iterator "const_iterator":
        #    PPL_Constraint operator*()
        #    PPL_const_iterator operator++()
        #    bint operator==(PPL_const_iterator&)
        #    bint operator!=(PPL_const_iterator&)
        #PPL_const_iterator begin()
        #PPL_const_iterator end()
        PPL_Constraint_System()
        PPL_Constraint_System(PPL_Constraint &g)
        PPL_Constraint_System(PPL_Constraint_System &gs)
        PPL_dimension_type space_dimension()
        bint has_equalities()
        bint has_strict_inequalities()
        void clear()
        void insert(PPL_Constraint &g)
        bint empty()
        void ascii_dump()
        bint OK()

    cdef enum PPL_Degenerate_Element "Parma_Polyhedra_Library::Degenerate_Element":
        UNIVERSE, EMPTY

    cdef enum PPL_Optimization_Mode "Parma_Polyhedra_Library::Optimization_Mode":
        MINIMIZATION, MAXIMIZATION

    cdef enum PPL_MIP_Problem_Status "Parma_Polyhedra_Library::MIP_Problem_Status":
        UNFEASIBLE_MIP_PROBLEM, UNBOUNDED_MIP_PROBLEM, OPTIMIZED_MIP_PROBLEM

    cdef cppclass PPL_Polyhedron:
        PPL_dimension_type space_dimension()
        PPL_dimension_type affine_dimension()
        PPL_Constraint_System& constraints()
        PPL_Constraint_System& minimized_constraints()
        PPL_Generator_System& generators()
        PPL_Generator_System& minimized_generators()
        PPL_Poly_Con_Relation relation_with(PPL_Constraint &c) except +ValueError
        PPL_Poly_Gen_Relation relation_with(PPL_Generator &g) except +ValueError
        bint is_empty()
        bint is_universe()
        bint is_topologically_closed()
        bint is_disjoint_from(PPL_Polyhedron &y) except +ValueError
        bint is_discrete()
        bint is_bounded()
        bint contains_integer_point()
        bint constrains(PPL_Variable var) except +ValueError
        bint bounds_from_above(PPL_Linear_Expression &expr) except +ValueError
        bint bounds_from_below(PPL_Linear_Expression &expr) except +ValueError
        bint maximize(PPL_Linear_Expression &expr, PPL_Coefficient &sup_n, PPL_Coefficient &sup_d,
                      cppbool &maximum)
        bint maximize(PPL_Linear_Expression &expr, PPL_Coefficient &sup_n, PPL_Coefficient &sup_d,
                      cppbool &maximum, PPL_Generator &g)
        bint minimize(PPL_Linear_Expression &expr, PPL_Coefficient &inf_n, PPL_Coefficient &inf_d,
                      cppbool &minimum)
        bint minimize(PPL_Linear_Expression &expr, PPL_Coefficient &inf_n, PPL_Coefficient &inf_d,
                      cppbool &minimum, PPL_Generator &g)
        bint frequency(PPL_Linear_Expression &expr, PPL_Coefficient &freq_n, PPL_Coefficient &freq_d,
                       PPL_Coefficient &val_n, PPL_Coefficient &val_d)
        bint contains(PPL_Polyhedron &y) except +ValueError
        bint strictly_contains(PPL_Polyhedron &y) except +ValueError
        void add_constraint(PPL_Constraint &c) except +ValueError
        void add_generator(PPL_Generator &g) except +ValueError
        void add_constraints(PPL_Constraint_System &cs) except +ValueError
        void add_generators(PPL_Generator_System &gs) except +ValueError
        void refine_with_constraint(PPL_Constraint &c) except +ValueError
        void refine_with_constraints(PPL_Constraint_System &cs) except +ValueError
        void unconstrain(PPL_Variable var) except +ValueError
        void intersection_assign(PPL_Polyhedron &y) except +ValueError
        void poly_hull_assign(PPL_Polyhedron &y) except +ValueError
        void upper_bound_assign(PPL_Polyhedron &y) except +ValueError
        void poly_difference_assign(PPL_Polyhedron &y) except +ValueError
        void difference_assign(PPL_Polyhedron &y) except +ValueError
        void drop_some_non_integer_points()
        void topological_closure_assign()
        void add_space_dimensions_and_embed(PPL_dimension_type m) except +ValueError
        void add_space_dimensions_and_project(PPL_dimension_type m) except +ValueError
        void concatenate_assign(PPL_Polyhedron &y) except +ValueError
        void remove_higher_space_dimensions(PPL_dimension_type new_dimension) except +ValueError
        void ascii_dump()
        int hash_code()
        PPL_dimension_type max_space_dimension()
        bint OK(cppbool check_not_empty=false)
        bint operator!=(PPL_Polyhedron &y)
        bint operator==(PPL_Polyhedron &y)

    cdef cppclass PPL_C_Polyhedron(PPL_Polyhedron):
        PPL_C_Polyhedron(PPL_dimension_type num_dimensions, PPL_Degenerate_Element)
        PPL_C_Polyhedron(PPL_Constraint_System &cs) except +ValueError
        PPL_C_Polyhedron(PPL_Generator_System &gs) except +ValueError
        PPL_C_Polyhedron(PPL_C_Polyhedron &y)

    cdef cppclass PPL_NNC_Polyhedron(PPL_Polyhedron):
        PPL_NNC_Polyhedron(PPL_dimension_type num_dimensions, PPL_Degenerate_Element kind)
        PPL_NNC_Polyhedron(PPL_Constraint_System &cs) except +ValueError
        PPL_NNC_Polyhedron(PPL_Generator_System &gs) except +ValueError
        PPL_NNC_Polyhedron(PPL_NNC_Polyhedron &y)
        PPL_NNC_Polyhedron(PPL_C_Polyhedron &y)

    cdef cppclass PPL_Poly_Gen_Relation:
        PPL_Poly_Gen_Relation(PPL_Poly_Gen_Relation &cpy_from)
        bint implies(PPL_Poly_Gen_Relation &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Poly_Con_Relation:
        PPL_Poly_Con_Relation(PPL_Poly_Con_Relation &cpy_from)
        bint implies(PPL_Poly_Con_Relation &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_MIP_Problem:
        PPL_MIP_Problem(PPL_MIP_Problem &cpy_from)
        PPL_MIP_Problem(PPL_dimension_type dim) except +ValueError
        PPL_MIP_Problem(PPL_dimension_type dim, PPL_Constraint_System &cs, PPL_Linear_Expression &obj, PPL_Optimization_Mode) except +ValueError
        PPL_dimension_type space_dimension()
        PPL_Linear_Expression& objective_function()
        void clear()
        void add_space_dimensions_and_embed(PPL_dimension_type m) except +ValueError
        void add_constraint(PPL_Constraint &c) except +ValueError
        void add_constraints(PPL_Constraint_System &cs) except +ValueError
        void set_objective_function(PPL_Linear_Expression &obj) except +ValueError
        void set_optimization_mode(PPL_Optimization_Mode mode)
        PPL_Optimization_Mode optimization_mode()
        bint is_satisfiable()
        PPL_MIP_Problem_Status solve()
        void evaluate_objective_function(PPL_Generator evaluating_point, PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        PPL_Generator& feasible_point()
        PPL_Generator optimizing_point() except +ValueError
        void optimal_value(PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        bint OK()
        PPL_MIP_Problem_Control_Parameter_Value get_control_parameter(PPL_MIP_Problem_Control_Parameter_Name name)
        void set_control_parameter(PPL_MIP_Problem_Control_Parameter_Value value)

cdef extern from "ppl.hh":
    PPL_Generator PPL_line          "Parma_Polyhedra_Library::line"          (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_ray           "Parma_Polyhedra_Library::ray"           (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_point         "Parma_Polyhedra_Library::point"         (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError
    PPL_Generator PPL_closure_point "Parma_Polyhedra_Library::closure_point" (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError


####################################################
# Cython does not support static methods; hack around
cdef extern from "ppl.hh":

    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_nothing  "Parma_Polyhedra_Library::Poly_Gen_Relation::nothing"  ()
    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_subsumes "Parma_Polyhedra_Library::Poly_Gen_Relation::subsumes" ()

    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_nothing  "Parma_Polyhedra_Library::Poly_Con_Relation::nothing"  ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_disjoint "Parma_Polyhedra_Library::Poly_Con_Relation::is_disjoint" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_strictly_intersects "Parma_Polyhedra_Library::Poly_Con_Relation::strictly_intersects" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_included "Parma_Polyhedra_Library::Poly_Con_Relation::is_included" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_saturates "Parma_Polyhedra_Library::Poly_Con_Relation::saturates" ()



####################################################
# Workaround for private constructors
cdef extern from "ppl_shim.hh":
    PPL_Poly_Gen_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Generator &g) except +ValueError
    PPL_Poly_Con_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Constraint &c) except +ValueError


####################################################
# C++ static methods not supported
# Note that the PPL_Generator default constructor is private, hence we must return pointers
cdef extern from "ppl_shim.hh":
    PPL_Generator* new_line(PPL_Linear_Expression &e) except +ValueError
    PPL_Generator* new_ray(PPL_Linear_Expression &e) except +ValueError
    PPL_Generator* new_point(PPL_Linear_Expression &e, PPL_Coefficient d) except +ValueError
    PPL_Generator* new_closure_point(PPL_Linear_Expression &e, PPL_Coefficient d) except +ValueError
    PPL_Generator* new_MIP_optimizing_point(PPL_MIP_Problem &problem) except +ValueError


cdef extern from "ppl_shim.hh":
    ctypedef void* gs_iterator_ptr
    cdef gs_iterator_ptr init_gs_iterator(PPL_Generator_System &gs)
    cdef PPL_Generator next_gs_iterator(gs_iterator_ptr)
    cdef bint is_end_gs_iterator(PPL_Generator_System &gs, gs_iterator_ptr gsi_ptr)
    cdef void delete_gs_iterator(gs_iterator_ptr)

    ctypedef void* cs_iterator_ptr
    cdef cs_iterator_ptr init_cs_iterator(PPL_Constraint_System &cs)
    cdef PPL_Constraint next_cs_iterator(cs_iterator_ptr)
    cdef bint is_end_cs_iterator(PPL_Constraint_System &cs, cs_iterator_ptr csi_ptr)
    cdef void delete_cs_iterator(cs_iterator_ptr)

    ctypedef void* mip_cs_iterator_ptr
    cdef mip_cs_iterator_ptr init_mip_cs_iterator(PPL_MIP_Problem& pb)
    cdef PPL_Constraint next_mip_cs_iterator(mip_cs_iterator_ptr mip_csi_ptr)
    cdef bint is_end_mip_cs_iterator(PPL_MIP_Problem &pb, mip_cs_iterator_ptr mip_csi_ptr)
    cdef void delete_mip_cs_iterator(mip_cs_iterator_ptr mip_csi_ptr)


cdef PPL_GeneratorType_str(PPL_GeneratorType t)

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

cdef class MIP_Problem(_mutable_or_immutable):
    cdef PPL_MIP_Problem *thisptr

cdef class MIP_Problem_constraints_iterator(object):
    cdef MIP_Problem pb
    cdef mip_cs_iterator_ptr mip_csi_ptr

cdef _wrap_Generator_System(PPL_Generator_System generator_system)
cdef _wrap_Constraint(PPL_Constraint constraint)
