from libcpp cimport bool as cppbool
from libcpp.vector cimport vector as cppvector

cdef extern from "gmp.h":
    # gmp integer
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

    void mpz_init(mpz_t)

cdef extern from "gmpxx.h":
    # gmp integer
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)

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

cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    ctypedef size_t PPL_dimension_type  "Parma_Polyhedra_Library::dimension_type"
    ctypedef mpz_class PPL_Coefficient  "Parma_Polyhedra_Library::Coefficient"
    cdef cppclass PPL_Variable          "Parma_Polyhedra_Library::Variable"
    cdef cppclass PPL_Variables_Set     "Parma_Polyhedra_Library::Variables_Set"
    cdef cppclass PPL_Linear_Expression "Parma_Polyhedra_Library::Linear_Expression"
    cdef cppclass PPL_Generator         "Parma_Polyhedra_Library::Generator"
    cdef cppclass PPL_Generator_System  "Parma_Polyhedra_Library::Generator_System"
    cdef cppclass PPL_Constraint        "Parma_Polyhedra_Library::Constraint"
    cdef cppclass PPL_Constraint_System "Parma_Polyhedra_Library::Constraint_System"
    cdef cppclass PPL_Congruence        "Parma_Polyhedra_Library::Congruence"
    cdef cppclass PPL_Congruence_System "Parma_Polyhedra_Library::Congruence_System"
    cdef cppclass PPL_Polyhedron        "Parma_Polyhedra_Library::Polyhedron"
    cdef cppclass PPL_C_Polyhedron      "Parma_Polyhedra_Library::C_Polyhedron"   (PPL_Polyhedron)
    cdef cppclass PPL_NNC_Polyhedron    "Parma_Polyhedra_Library::NNC_Polyhedron" (PPL_Polyhedron)
    cdef cppclass PPL_Poly_Gen_Relation "Parma_Polyhedra_Library::Poly_Gen_Relation"
    cdef cppclass PPL_Poly_Con_Relation "Parma_Polyhedra_Library::Poly_Con_Relation"
    cdef cppclass PPL_MIP_Problem       "Parma_Polyhedra_Library::MIP_Problem"
    cdef cppclass PPL_mip_iterator      "Parma_Polyhedra_Library::MIP_Problem::const_iterator"
    cdef cppclass PPL_gs_iterator       "Parma_Polyhedra_Library::Generator_System::const_iterator"
    cdef cppclass PPL_Constraint_System_iterator   "Parma_Polyhedra_Library::Constraint_System::const_iterator"
    cdef cppclass PPL_Congruence_System_iterator   "Parma_Polyhedra_Library::Congruence_System::const_iterator"

    cdef cppclass PPL_Bit_Row           "Parma_Polyhedra_Library::Bit_Row"
    cdef cppclass PPL_Bit_Matrix        "Parma_Polyhedra_Library::Bit_Matrix"

    cdef cppclass PPL_Variable:
        PPL_Variable(PPL_dimension_type i)
        PPL_dimension_type id()
        PPL_dimension_type space_dimension()

    cdef cppclass PPL_Variables_Set:
        PPL_Variables_Set()
        PPL_Variables_Set(PPL_Variable v)
        PPL_Variables_Set(PPL_Variable v, PPL_Variable w)
        PPL_dimension_type space_dimension()
        void insert(PPL_Variable v)
        size_t size()
        void ascii_dump()

    # class Parma_Polyhedra_Library::Linear_Expression
    # lines 28238-28879 of ppl.hh
    cdef cppclass PPL_Linear_Expression:
        PPL_Linear_Expression()
        PPL_Linear_Expression(PPL_Linear_Expression &e)
        PPL_Linear_Expression(PPL_Coefficient n)
        PPL_Linear_Expression(PPL_Variable v)

        PPL_dimension_type max_space_dimension()
        PPL_dimension_type space_dimension()
        void set_space_dimension(PPL_dimension_type d)
        PPL_Coefficient coefficient(PPL_Variable v)
        void set_coefficient(PPL_Variable v, PPL_Coefficient)
        PPL_Coefficient inhomogeneous_term()
        void set_inhomogeneous_term(PPL_Coefficient n)
        void linear_combine(const PPL_Linear_Expression& y, PPL_Variable v)
        void linear_combine(const PPL_Linear_Expression& y, PPL_Coefficient c1, PPL_Coefficient c2)
        void linear_combine_lax(const PPL_Linear_Expression& u, PPL_Coefficient c1, PPL_Coefficient c2)
        void swap_space_dimensions(PPL_Variable v1, PPL_Variable v2)
        void remove_space_dimensions(const PPL_Variables_Set)
        void shift_space_dimensions(PPL_Variable v, PPL_dimension_type n)
        void permute_space_dimensions(const cppvector[PPL_Variable]& cycle) except +ValueError
        bint is_zero()
        bint all_homogeneous_terms_are_zero()
        bint is_equal_to(PPL_Linear_Expression& x)
        bint all_zeroes(const PPL_Variables_Set& v)

        void ascii_dump()

        #PPL_Linear_Expression operator+=(PPL_Linear_Expression& e)
        #PPL_Linear_Expression operator-=(PPL_Linear_Expression& e)
        #PPL_Linear_Expression operator*=(PPL_Coefficient n)
        #PPL_Linear_Expression operator/=(PPL_Coefficient n)
        PPL_Linear_Expression operator+(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator-(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator*(PPL_Coefficient n)
        PPL_Constraint operator> (PPL_Linear_Expression& e)
        PPL_Constraint operator>=(PPL_Linear_Expression& e)
        PPL_Constraint operator==(PPL_Linear_Expression& e)
        PPL_Constraint operator<=(PPL_Linear_Expression& e)
        PPL_Constraint operator< (PPL_Linear_Expression& e)

    cdef cppclass PPL_Constraint:
        PPL_Constraint()
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
        void permute_space_dimensions(const cppvector[PPL_Variable]& cycle) except +ValueError

    cdef cppclass PPL_Generator:
        PPL_Generator(PPL_Generator &g)
        PPL_dimension_type space_dimension()
        void set_space_dimension(PPL_dimension_type n)
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
        void permute_space_dimensions(const cppvector[PPL_Variable]& cycle) except +ValueError

    cdef cppclass PPL_Congruence:
        PPL_Congruence()
        PPL_Congruence(const PPL_Congruence &g)
        PPL_Congruence(const PPL_Constraint &c) except +ValueError
        PPL_dimension_type space_dimension()
        # NOTE: curiously, this can raise an error (behavior different from Linear_Expression)
        PPL_Coefficient coefficient(PPL_Variable v) except +ValueError
        PPL_Coefficient inhomogeneous_term()
        PPL_Coefficient modulus()
        void set_modulus(PPL_Coefficient& m)
        void scale(PPL_Coefficient& m)
        bint is_tautological()
        bint is_inconsistent()
        bint is_proper_congruence()
        bint is_equality()
        void ascii_dump()
        void swap_space_dimension(PPL_Variable v1, PPL_Variable v2)
        void set_space_dimension(PPL_dimension_type n)
        void shift_space_dimensions(PPL_Variable v, PPL_dimension_type n)
        void sign_normalize()
        void normalize()
        void strong_normalize()
        PPL_dimension_type max_space_dimension()

        cppbool operator==(const PPL_Congruence &x, const PPL_Congruence &y)
        cppbool operator!=(const PPL_Congruence &x, const PPL_Congruence &y)

    cdef cppclass PPL_Congruence_System:
        PPL_Congruence_System()
        PPL_Congruence_System(PPL_Congruence &c)
        PPL_Congruence_System(PPL_Congruence_System &cs)
        PPL_dimension_type space_dimension()
        PPL_Congruence_System_iterator begin()
        PPL_Congruence_System_iterator end()
        bint has_equalities()
        bint has_strict_inequalities()
        void clear()
        void insert(PPL_Congruence &g)
        bint empty()
        void ascii_dump()

    cdef cppclass PPL_Congruence_System_iterator:
        PPL_Congruence_System_iterator()
        PPL_Congruence_System_iterator(PPL_Congruence_System_iterator &csi)
        PPL_Congruence& operator* ()
        PPL_Congruence_System_iterator inc "operator++" (int i)
        cppbool operator==(PPL_Congruence_System_iterator& y)
        cppbool operator!=(PPL_Congruence_System_iterator& y)

    cdef cppclass PPL_Generator_System:
        PPL_Generator_System()
        PPL_Generator_System(PPL_Generator &g)
        PPL_Generator_System(PPL_Generator_System &gs)
        PPL_dimension_type space_dimension()
        void set_space_dimension(PPL_dimension_type space_dim)
        PPL_gs_iterator begin()
        PPL_gs_iterator end()
        void clear()
        void insert(PPL_Generator &g)
        bint empty()
        void ascii_dump()

    cdef cppclass PPL_mip_iterator:
        PPL_mip_iterator(PPL_mip_iterator &mipi)
        PPL_Constraint& operator* ()
        PPL_mip_iterator inc "operator++" (int i)
        cppbool operator==(PPL_mip_iterator& y)
        cppbool operator!=(PPL_mip_iterator& y)

    cdef cppclass PPL_gs_iterator:
        PPL_gs_iterator()
        PPL_gs_iterator(PPL_gs_iterator &gsi)
        PPL_Generator& operator* ()
        PPL_gs_iterator inc "operator++" (int i)
        cppbool operator==(PPL_gs_iterator& y)
        cppbool operator!=(PPL_gs_iterator& y)

    cdef cppclass PPL_Constraint_System_iterator:
        PPL_Constraint_System_iterator()
        PPL_Constraint_System_iterator(PPL_Constraint_System_iterator &csi)
        PPL_Constraint& operator* ()
        PPL_Constraint_System_iterator inc "operator++" (int i)
        cppbool operator==(PPL_Constraint_System_iterator& y)
        cppbool operator!=(PPL_Constraint_System_iterator& y)

    cdef cppclass PPL_Constraint_System:
        PPL_Constraint_System()
        PPL_Constraint_System(PPL_Constraint &g)
        PPL_Constraint_System(PPL_Constraint_System &gs)
        PPL_dimension_type space_dimension()
        PPL_Constraint_System_iterator begin()
        PPL_Constraint_System_iterator end()
        bint has_equalities()
        bint has_strict_inequalities()
        void clear()
        void insert(PPL_Constraint &g)
        bint empty()
        void ascii_dump()

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
        void BHRZ03_widening_assign(PPL_Polyhedron &y, unsigned* tp) except +ValueError
        void limited_BHRZ03_extrapolation_assign(PPL_Polyhedron &y, PPL_Constraint_System &cs, unsigned* tp) except +ValueError
        void bounded_BHRZ03_extrapolation_assign(PPL_Polyhedron &y, PPL_Constraint_System &cs, unsigned* tp) except +ValueError
        void H79_widening_assign(PPL_Polyhedron &y, unsigned* tp) except +ValueError
        void widening_assign(PPL_Polyhedron &y, unsigned* tp) except +ValueError
        void limited_H79_extrapolation_assign(PPL_Polyhedron &y, PPL_Constraint_System &cs, unsigned* tp) except +ValueError
        void bounded_H79_extrapolation_assign(PPL_Polyhedron &y, PPL_Constraint_System &cs, unsigned* tp) except +ValueError
        void add_space_dimensions_and_embed(PPL_dimension_type m) except +ValueError
        void add_space_dimensions_and_project(PPL_dimension_type m) except +ValueError
        void concatenate_assign(PPL_Polyhedron &y) except +ValueError
        void remove_higher_space_dimensions(PPL_dimension_type new_dimension) except +ValueError
        void affine_image(const PPL_Variable, const PPL_Linear_Expression& expr) except +ValueError
        void affine_preimage(const PPL_Variable, const PPL_Linear_Expression& expr) except +ValueError
        void ascii_dump()
        int hash_code()
        PPL_dimension_type max_space_dimension()
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

    cdef cppclass PPL_Poly_Con_Relation:
        PPL_Poly_Con_Relation(PPL_Poly_Con_Relation &cpy_from)
        bint implies(PPL_Poly_Con_Relation &y)
        void ascii_dump()

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
        void add_to_integer_space_dimensions(PPL_Variables_Set &i_vars) except +ValueError
        void set_objective_function(PPL_Linear_Expression &obj) except +ValueError
        void set_optimization_mode(PPL_Optimization_Mode mode)
        PPL_Optimization_Mode optimization_mode()
        bint is_satisfiable()
        PPL_MIP_Problem_Status solve()
        void evaluate_objective_function(PPL_Generator evaluating_point, PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        PPL_Generator& feasible_point()
        PPL_Generator optimizing_point() except +ValueError
        void optimal_value(PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        PPL_MIP_Problem_Control_Parameter_Value get_control_parameter(PPL_MIP_Problem_Control_Parameter_Name name)
        void set_control_parameter(PPL_MIP_Problem_Control_Parameter_Value value)
        PPL_mip_iterator constraints_begin()
        PPL_mip_iterator constraints_end()

    cdef cppclass PPL_Bit_Row:
        PPL_Bit_Row()
        PPL_Bit_Row(const PPL_Bit_Row& y, const PPL_Bit_Row& z)
        void set(unsigned long k)
        void set_until(unsigned long k)
        void clear_from(unsigned long k)
        void clear()
        void union_assign(const PPL_Bit_Row& x, const PPL_Bit_Row& y)
        void intersection_assign(const PPL_Bit_Row& x, const PPL_Bit_Row& y)
        void difference_assign(const PPL_Bit_Row&x, const PPL_Bit_Row& y)

        unsigned long first()
        unsigned long last()
        unsigned long prev(unsigned long position)
        unsigned long next(unsigned long position)
        unsigned long count_ones()
        cppbool empty()

    cdef cppclass PPL_Bit_Matrix:
        PPL_Bit_Matrix()
        PPL_Bit_Matrix(PPL_dimension_type n_rows, PPL_dimension_type n_columns)
        PPL_Bit_Matrix(const PPL_Bit_Matrix& y)

        PPL_Bit_Row& operator[](PPL_dimension_type k)
        const PPL_Bit_Row& operator[](PPL_dimension_type k)

        void transpose()
        void transpose_assign(const PPL_Bit_Matrix& y)

        PPL_dimension_type num_columns()
        PPL_dimension_type num_rows()

        void sort_rows()

cdef extern from "ppl.hh":
    PPL_Generator PPL_line          "Parma_Polyhedra_Library::line"             (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_ray           "Parma_Polyhedra_Library::ray"              (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_point         "Parma_Polyhedra_Library::point"            (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError
    PPL_Generator PPL_closure_point "Parma_Polyhedra_Library::closure_point"    (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError

cdef extern from "ppl.hh":

    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_nothing  "Parma_Polyhedra_Library::Poly_Gen_Relation::nothing"  ()
    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_subsumes "Parma_Polyhedra_Library::Poly_Gen_Relation::subsumes" ()

    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_nothing  "Parma_Polyhedra_Library::Poly_Con_Relation::nothing"  ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_disjoint "Parma_Polyhedra_Library::Poly_Con_Relation::is_disjoint" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_strictly_intersects "Parma_Polyhedra_Library::Poly_Con_Relation::strictly_intersects" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_included "Parma_Polyhedra_Library::Poly_Con_Relation::is_included" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_saturates "Parma_Polyhedra_Library::Poly_Con_Relation::saturates" ()

cdef extern from "ppl_shim.hh":
    PPL_Poly_Gen_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Generator &g) except +ValueError
    PPL_Poly_Con_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Constraint &c) except +ValueError

    PPL_Congruence modulo(const PPL_Linear_Expression &expr, PPL_Coefficient& mod)
