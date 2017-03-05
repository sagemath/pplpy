import ppl
import doctest
import sys
res = doctest.testmod(ppl, optionflags=doctest.ELLIPSIS)
print(res)
sys.exit(res[0])