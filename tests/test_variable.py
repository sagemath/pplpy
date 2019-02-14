import unittest

from ppl import Variable

class TestVariable(unittest.TestCase):
    def test_creation_valid(self):
        Variable(0)
    
    def test_creation_invalid(self):
        self.assertRaises(OverflowError, Variable, -1)
        self.assertRaises(TypeError, Variable, "hello")
