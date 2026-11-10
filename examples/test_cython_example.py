from gmpy2 import mpq

from ppl import C_Polyhedron, Generator_System, Variable, point, ray

import cython_example


def test_mip_problem():
    m = cython_example.mip_problem()
    assert m.optimal_value() == mpq(10, 3)
    p = m.optimizing_point()
    assert p.coefficients() == (10, 0)
    assert p.divisor() == 3


def test_polyhedron():
    x, y, z = Variable(0), Variable(1), Variable(2)
    gs = Generator_System()
    gs.insert(point(z))
    gs.insert(ray(x - z))
    gs.insert(ray(y - z))
    assert cython_example.polyhedron() == C_Polyhedron(gs)


def test_space_dimension():
    assert cython_example.space_dimension(cython_example.polyhedron()) == 3
