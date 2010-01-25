#!/usr/bin/env python
"""Somewhat hacky unittest runner that assumes all scripts in the unittests directory are unittests
and have a class named Test containing all their tests. It then runs them all as a single suite.
"""

import os
import unittest
import unittests
from unittests import *

# some unittests expect files existing in the unittests dir
os.chdir('unittests')

suite = unittest.TestSuite()
for test in unittests.__all__:
    clss = eval(test+".Test")
    s = unittest.TestLoader().loadTestsFromTestCase(clss)
    suite.addTest(s)

unittest.TextTestRunner(verbosity=2).run(suite)
