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

from ppl.wrappers import _mutable_or_immutable as ExampleObj
from ppl.wrappers cimport Variable, Constraint_System, MIP_Problem
"""
By default, any object is mutable. It can then be `set_immutable`.
"""
def f():
    x = ExampleObj()
    assert(x.is_mutable())
    assert(not x.is_immutable())
    x.set_immutable()
    assert(x.is_immutable())

    x = Variable(0)
    y = Variable(1)
    cs = Constraint_System()
    cs.insert( x >= 0)
    cs.insert( y >= 0 )
    cs.insert( 3 * x + 5 * y <= 10 )
    m = MIP_Problem(2, cs, x + y)
    m.objective_function()
    #x0+x1
    print("cython call of pplpy is a success")
f()