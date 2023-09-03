# distutils: language = c++
# distutils: libraries = gmp gmpxx ppl m
#*****************************************************************************
#       Copyright (C) 2010-2014 Volker Braun  <vbraun.name@gmail.com>
#                     2011 Simon King <simon.king@uni-jena.de>
#                     2011-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2012 Risan <ptrrsn.1@gmail.com>
#                     2013 Julien Puydt <julien.puydt@laposte.net>
#                     2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#                     2015 André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                     2016 Jori Mäntysalo <jori.mantysalo@uta.fi>
#                     2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#                     2016-2017 Frédéric Chapoton <chapoton@math.univ-lyon1.fr>
#                     2016-2018 Vincent Delecroix <vincent.delecroix@labri.fr>
#                     2017-2018 Vincent Klein <vinklein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from gmpy2 cimport import_gmpy2, GMPy_MPQ_From_mpz
from cython.operator cimport dereference as deref

# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# initialize gmpy2 C API
import_gmpy2()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()


cdef class MIP_Problem(object):
    r"""
    wrapper for PPL's MIP_Problem class

    An object of the class MIP_Problem represents a Mixed Integer
    (Linear) Program problem.

    INPUT:

    - ``dim`` -- integer
    - ``args`` -- an array of the defining data of the MIP_Problem.
      For each element, any one of the following is accepted:

      * A :class:`Constraint_System`.

      * A :class:`Linear_Expression`.

    OUTPUT:

    A :class:`MIP_Problem`.

    Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3 * x + 5 * y <= 10)
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m
        MIP Problem (maximization, 2 variables, 3 constraints)
        >>> m.optimal_value()
        mpq(10,3)
        >>> float(_)
        3.333333333333333
        >>> m.optimizing_point()
        point(10/3, 0/3)
    """
    def __repr__(self):
        """
        String representation of MIP Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m
        MIP Problem (maximization, 2 variables, 3 constraints)
        """
        dim = self.space_dimension()
        ncs = sum(1 for _ in self)
        return 'MIP Problem ({}, {} variable{}, {} constraint{})'.format(
                self.optimization_mode(),
                dim,
                's' if dim > 1 else '',
                ncs,
                's' if ncs > 1 else '')

    def __cinit__(self, PPL_dimension_type dim = 0, *args):
        """
        The Cython constructor.

        Tests:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> MIP_Problem(0)
        A MIP_Problem
        Maximize: 0
        Subject to constraints
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x + y <= 2)
        >>> M = MIP_Problem(2, cs, 0)
        >>> M = MIP_Problem(2, cs, x)
        >>> M = MIP_Problem(2, None, None)
        Traceback (most recent call last):
        ...
        TypeError: Cannot convert NoneType to ppl.Constraint_System
        >>> M = MIP_Problem(2, cs, 'hey')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'hey' to an integer
        >>> M = MIP_Problem(2, cs, x, 'middle')
        Traceback (most recent call last):
        ...
        ValueError: unknown mode 'middle'
        """
        if not args:
            self.thisptr = new PPL_MIP_Problem(dim)
            return
        elif len(args) == 1:
            raise ValueError('cannot initialize from {}'.format(args))

        cdef Constraint_System cs = <Constraint_System?>args[0]
        cdef Linear_Expression obj
        try:
            obj = <Linear_Expression?> args[1]
        except TypeError:
            obj = Linear_Expression(args[1])

        cdef PPL_Optimization_Mode mode = MAXIMIZATION
        if len(args) == 3:
            if args[2] == 'maximization':
                mode = MAXIMIZATION
            elif args[2] == 'minimization':
                mode = MINIMIZATION
            else:
                raise ValueError('unknown mode {!r}'.format(args[2]))
        self.thisptr = new PPL_MIP_Problem(dim, cs.thisptr[0], obj.thisptr[0], mode)

    def __dealloc__(self):
        """
        The Cython destructor
        """
        del self.thisptr

    def __iter__(self):
        r"""
        Iterator through the constraints
        Tests:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> M = MIP_Problem(2)
        >>> for c in M: print(c)
        >>> M.add_constraint(x + y <= 5)
        >>> for c in M: print(c)
        -x0-x1+5>=0
        >>> M.add_constraint(3*x - 18*y >= -2)
        >>> for c in M: print(c)
        -x0-x1+5>=0
        3*x0-18*x1+2>=0
        >>> M = MIP_Problem(1)
        >>> M.add_constraint(x <= 5)
        >>> it = M.__iter__()
        >>> next(it)
        -x0+5>=0
        >>> next(it)
        Traceback (most recent call last):
        ...
        StopIteration
        """
        cdef PPL_mip_iterator *mip_csi_ptr = new PPL_mip_iterator(self.thisptr[0].constraints_begin())
        try:
            while mip_csi_ptr[0] != self.thisptr[0].constraints_end():
                yield _wrap_Constraint(deref(mip_csi_ptr[0].inc(1)))
        finally:
            del mip_csi_ptr

    def constraints(self):
        r"""
        Return the constraints of this MIP

        The output is an instance of :class:`Constraint_System`.

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> M = MIP_Problem(2)
        >>> M.add_constraint(x + y <= 5)
        >>> M.add_constraint(3*x - 18*y >= -2)
        >>> M.constraints()
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0}

        Note that modifying the output of this method will not modify the
        underlying MIP problem object:

        >>> cs = M.constraints()
        >>> cs.insert(x <= 3)
        >>> cs
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0, -x0+3>=0}
        >>> M.constraints()
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0}
        """
        cdef Constraint_System c = Constraint_System(None)
        cdef PPL_Constraint_System* cs = new PPL_Constraint_System()
        cdef PPL_mip_iterator* mip_it = new PPL_mip_iterator(self.thisptr[0].constraints_begin())

        while mip_it[0] != self.thisptr[0].constraints_end():
            cs[0].insert(deref(mip_it[0]))
            mip_it[0].inc(1)
        c.thisptr = cs
        del mip_it
        return c

    def optimization_mode(self):
        """
        Return the optimization mode used in the MIP_Problem.

        It will return "maximization" if the MIP_Problem was set
        to MAXIMIZATION mode, and "minimization" otherwise.

        Examples:

            >>> from ppl import MIP_Problem
            >>> m = MIP_Problem()
            >>> m.optimization_mode()
            'maximization'
        """
        if self.thisptr.optimization_mode() == MAXIMIZATION:
            return "maximization"
        elif self.thisptr.optimization_mode() == MINIMIZATION:
            return "minimization"

    def optimal_value(self):
        """
        Return the optimal value of the MIP_Problem. ValueError thrown if self does not
        have an optimizing point, i.e., if the MIP problem is unbounded or not satisfiable.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.optimal_value()
        mpq(10,3)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> m = MIP_Problem(1, cs, x + x )
        >>> m.optimal_value()
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::optimizing_point():
        *this ... have an optimizing point.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.optimal_value(sup_n, sup_d)
        finally:
            sig_off()
        return GMPy_MPQ_From_mpz(sup_n.get_mpz_t(), sup_d.get_mpz_t())

    def space_dimension(self):
        """
        Return the space dimension of the MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def objective_function(self):
        """
        Return the optimal value of the MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.objective_function()
        x0+x1
        """
        rc = Linear_Expression()
        rc.thisptr[0] = self.thisptr.objective_function()
        return rc

    def clear(self):
        """
        Reset the MIP_Problem to be equal to the trivial MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.objective_function()
        x0+x1
        >>> m.clear()
        >>> m.objective_function()
        0
        """
        self.thisptr.clear()

    def add_space_dimensions_and_embed(self, PPL_dimension_type m):
        """
        Adds m new space dimensions and embeds the old MIP problem in the new vector space.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.add_space_dimensions_and_embed(5)
        >>> m.space_dimension()
        7
        """
        sig_on()
        self.thisptr.add_space_dimensions_and_embed(m)
        sig_off()

    def add_constraint(self, Constraint c):
        """
        Adds a copy of constraint c to the MIP problem.

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimal_value()
        mpq(10,3)

        Tests:

        >>> z = Variable(2)
        >>> m.add_constraint(z >= -3)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::add_constraint(c):
        c.space_dimension() == 3 exceeds this->space_dimension == 2.
        """
        sig_on()
        try:
            self.thisptr.add_constraint(c.thisptr[0])
        finally:
            sig_off()

    def add_constraints(self, Constraint_System cs):
        """
        Adds a copy of the constraints in cs to the MIP problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3 * x + 5 * y <= 10)
        >>> m = MIP_Problem(2)
        >>> m.set_objective_function(x + y)
        >>> m.add_constraints(cs)
        >>> m.optimal_value()
        mpq(10,3)

        Tests:

        >>> p = Variable(9)
        >>> cs.insert(p >= -3)
        >>> m.add_constraints(cs)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::add_constraints(cs):
        cs.space_dimension() == 10 exceeds this->space_dimension() == 2.
        """
        sig_on()
        try:
            self.thisptr.add_constraints(cs.thisptr[0])
        finally:
            sig_off()

    def add_to_integer_space_dimensions(self, Variables_Set i_vars):
        """
        Sets the variables whose indexes are in set `i_vars` to be integer space dimensions.

        Examples:

        >>> from ppl import Variable, Variables_Set, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2)
        >>> m.set_objective_function(x + y)
        >>> m.add_constraints(cs)
        >>> i_vars = Variables_Set(x, y)
        >>> m.add_to_integer_space_dimensions(i_vars)
        >>> m.optimal_value()
        mpq(3,1)
        """
        sig_on()
        try:
            self.thisptr.add_to_integer_space_dimensions(i_vars.thisptr[0])
        finally:
            sig_off()

    def set_objective_function(self, obj):
        """
        Sets the objective function to obj.

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimal_value()
        mpq(10,3)

        Tests:

        >>> z = Variable(2)
        >>> m.set_objective_function(x + y + z)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::set_objective_function(obj):
        obj.space_dimension() == 3 exceeds this->space_dimension == 2.

        >>> M = MIP_Problem(1)
        >>> M.set_objective_function(Variable(0))
        """
        if isinstance(obj, Variable):
            obj = Linear_Expression(obj)
        elif not isinstance(obj, Linear_Expression):
            raise ValueError('not an objective function')
        self.thisptr.set_objective_function((<Linear_Expression> obj).thisptr[0])

    def set_optimization_mode(self, mode):
        """
        Sets the optimization mode to mode.

        Examples:

        >>> from ppl import MIP_Problem
        >>> m = MIP_Problem()
        >>> m.optimization_mode()
        'maximization'
        >>> m.set_optimization_mode('minimization')
        >>> m.optimization_mode()
        'minimization'

        Tests:

        >>> m.set_optimization_mode('max')
        Traceback (most recent call last):
        ...
        ValueError: Unknown value: mode=max.
        """
        if mode == 'minimization':
            self.thisptr.set_optimization_mode(MINIMIZATION)
        elif mode == 'maximization':
            self.thisptr.set_optimization_mode(MAXIMIZATION)
        else:
            raise ValueError('Unknown value: mode='+str(mode)+'.')

    def is_satisfiable(self):
        """
        Check if the MIP_Problem is satisfiable

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.is_satisfiable()
        True
        """
        return self.thisptr.is_satisfiable()

    def evaluate_objective_function(self, Generator evaluating_point):
        """
        Return the result of evaluating the objective function on evaluating_point. ValueError thrown
        if self and evaluating_point are dimension-incompatible or if the generator evaluating_point is not a point.

        Examples:

        >>> from ppl import Variable, MIP_Problem, Generator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> g = Generator.point(5 * x - 2 * y, 7)
        >>> m.evaluate_objective_function(g)
        mpq(3,7)
        >>> z = Variable(2)
        >>> g = Generator.point(5 * x - 2 * z, 7)
        >>> m.evaluate_objective_function(g)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::evaluate_objective_function(p, n, d):
        *this and p are dimension incompatible.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.evaluate_objective_function(evaluating_point.thisptr[0], sup_n, sup_d)
        finally:
            sig_off()

        return GMPy_MPQ_From_mpz(sup_n.get_mpz_t(), sup_d.get_mpz_t())

    def solve(self):
        """
        Optimizes the MIP_Problem

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.solve()
        {'status': 'optimized'}
        """
        sig_on()
        try:
            tmp = self.thisptr.solve()
        finally:
            sig_off()
        if tmp == UNFEASIBLE_MIP_PROBLEM:
            return {'status': 'unfeasible'}
        elif tmp == UNBOUNDED_MIP_PROBLEM:
            return {'status': 'unbounded'}
        else:
            return {'status': 'optimized'}

    def optimizing_point(self):
        """
        Returns an optimal point for the MIP_Problem, if it exists.

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimizing_point()
        point(10/3, 0/3)
        """
        cdef PPL_Generator *g
        sig_on()
        try:
            g = new PPL_Generator(self.thisptr[0].optimizing_point())
        finally:
            sig_off()
        return _wrap_Generator(g[0])
