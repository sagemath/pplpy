# distutils: language = c++
# distutils: libraries = ppl
"""
The goal of this file is to test cython can use pplpy package properly
In order to do this we do some test with objects from each packages and extension :
ppl
ppl.linear_algebra
ppl.constraint
ppl.mip_problem
"""

from ppl.linear_algebra cimport Variable
from ppl.constraint cimport Constraint_System
from ppl.mip_problem cimport MIP_Problem
from ppl.polyhedron cimport C_Polyhedron

def test():
    x = Variable(0)
    y = Variable(1)
    cs = Constraint_System()
    cs.insert( x >= 0)
    cs.insert( y >= 0 )
    cs.insert( 3 * x + 5 * y <= 10 )
    m = MIP_Problem(2, cs, x + y)
    m.objective_function()

    from ppl import C_Polyhedron
    C_Polyhedron( 5*x-2*y >=  x+y-1 )

    print("-"*80)
    print("Cython test 1 OK")
    print("-"*80)

def example():
    "Cython version of the example from the README"
    cdef Variable x = Variable(0)
    cdef Variable y = Variable(1)
    cdef Variable z = Variable(2)
    cdef Constraint_System cs = Constraint_System()
    cs.insert(x >= 0)
    cs.insert(y >= 0)
    cs.insert(x + y + z == 1)
    cdef C_Polyhedron poly = C_Polyhedron(cs)
    print(poly.minimized_generators())
    print('dim = %lu' % poly.thisptr.space_dimension())
    print("-"*80)
    print("Cython example 1 OK")
    print("-"*80)
