import ppl
import doctest
import sys
import subprocess

#  - python setup.py build_ext --inplace
#  - python -c "import testpplpy"
#  - python setup2.py build_ext --inplace
#  - python -c "import testpplpy2"
subprocess.call([sys.executable, 'setup.py', 'build_ext', '--inplace'])
subprocess.call([sys.executable, '-c', '"import testpplpy"'])
subprocess.call([sys.executable, 'setup2.py', 'build_ext', '--inplace'])
subprocess.call([sys.executable, '-c', '"import testpplpy2"'])

print("Running doctests for pplpy")
print('-'*80)
ans = 0
for mod in [ppl, ppl.wrappers]:
    res = doctest.testmod(mod, optionflags=doctest.ELLIPSIS | doctest.REPORT_NDIFF)
    print(mod)
    print(res)
    print('-'*80)
    ans = ans | res[0]
sys.exit(ans)
