import unittest

from ppl import Variable, Constraint, Constraint_System

class TestConstraint(unittest.TestCase):
    def test_creation_empty(self):
        c = Constraint()
        self.assertTrue(c.is_tautological())

    def test_creation_other(self):
        x = Variable(0)
        y = Variable(1)
        c = x + 3 * y == 1
        cc = Constraint(c)
        self.assertTrue(c.is_equivalent_to(cc))

    def test_creation_invalid(self):
        self.assertRaises(TypeError, Constraint, "hello")

class TestConstraint_System(unittest.TestCase):
    def test_creation_empty(self):
        cs = Constraint_System()
        self.assertTrue(cs.empty())

    def test_creation_other(self):
        x = Variable(0)
        y = Variable(1)
        cs = Constraint_System(5*x - 2*y > 0)
        cs.insert(6 * x < 3 * y)
        ccs = Constraint_System(cs)
        self.assertTrue(cs[0].is_equivalent_to(ccs[0]))

    def test_creation_invalid(self):
        self.assertRaises(TypeError, Constraint_System, "hello")
