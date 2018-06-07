import ppl
import doctest
import sys
import os

path = os.path.dirname(__file__)
if path:
    os.chdir(path)

print("Running doctests for pplpy")
print('-'*80)
ans = 0
for mod in [ppl, ppl.linear_algebra, ppl.mip_problem, ppl.polyhedron, ppl.generator, ppl.constraint, ppl.bit_arrays]:
    res = doctest.testmod(mod, optionflags=doctest.ELLIPSIS | doctest.REPORT_NDIFF)
    print(mod)
    print(res)
    print('-'*80)
    ans = ans | res[0]
sys.exit(ans)
