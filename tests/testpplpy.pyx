# distutils: language = c++
# distutils: libraries = ppl
"""
The goal of this file is to test cython can use pplpy package properly
In order to do this we do some test with objects from each packages and extension :
ppl
ppl.wrappers
ppl.cygmp
ppl.cygmp.pylong
ppl.cygmp.utils
"""

from ppl.wrappers cimport Variable, Constraint_System
from ppl.mip_problem cimport MIP_Problem

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
