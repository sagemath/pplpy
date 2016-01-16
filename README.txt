==================
PPL Python wrapper
==================

This Python package provides a wrapper to the C++ [Parma Polyhedra Library
(PPL)](http://bugseng.com/products/ppl/). The names of objects and methods
are the same as in the library:

    import ppl
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    z = ppl.Variable(2)
    cs = Constraint_System()
    cs.insert(x >= 0)
    cs.insert(y >= 0)
    cs.insert(z >= 0)
    cs.insert(x + y + z == 1)
    poly = C_Polyhedron(cs)

Note that in order to compile this package you need the ppl libraries and
header files. 
