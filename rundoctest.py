import ppl
import doctest
import sys
res = str(doctest.testmod(ppl, optionflags=doctest.ELLIPSIS))
print(res)
offset = res.find("TestResults(failed=") + len("TestResults(failed=")
sys.exit(int(res[offset:res.find(',', offset)]))
