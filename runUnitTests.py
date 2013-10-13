#!/usr/bin/env python
"""Somewhat hacky unittest runner that assumes all scripts in the unittests directory are unittests
and have a class containing all their tests. It then runs them all as a single suite.
"""

import os
import inspect
import unittest
import unittests
from unittests import *

suite = unittest.TestSuite()
for mod in inspect.getmembers(unittests, inspect.ismodule):
    s = unittest.TestLoader().loadTestsFromModule(mod[1])
    suite.addTest(s)

# some unittests expect files existing in the unittests dir
os.chdir('unittests')
unittest.TextTestRunner(verbosity=2).run(suite)
