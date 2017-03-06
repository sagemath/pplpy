import ppl
import doctest
import sys
res = doctest.testmod(ppl.ppl, optionflags=doctest.ELLIPSIS|doctest.REPORT_NDIFF)
print(res)
sys.exit(res[0])