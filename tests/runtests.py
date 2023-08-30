import ppl
import sys
import os

path = os.path.dirname(__file__)
if path:
    os.chdir(path)

ans = 0 # set to nonzero if an error is found

print("Running pplpy doctests")
print('-'*80)
import doctest
for mod in [ppl, ppl.linear_algebra, ppl.mip_problem, ppl.polyhedron, ppl.generator, ppl.constraint, ppl.congruence, ppl.bit_arrays]:
    res = doctest.testmod(mod, optionflags=doctest.ELLIPSIS | doctest.REPORT_NDIFF | doctest.NORMALIZE_WHITESPACE)
    print(mod)
    print(res)
    print('-'*80)
    ans = ans | res[0]

print("Running unittests")
print('-'*80)
import unittest
for mod in ['test_variable', 'test_constraint']:
    res = unittest.main(module=mod, exit=False, failfast=False, verbosity=2)
    ans = ans | bool(res.result.errors)

sys.exit(ans)
