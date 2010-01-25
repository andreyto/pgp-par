"Unit test package, auto import all .py programs as unit tests"
import os

__all__ = [os.path.splitext(x)[0] for x in os.listdir('unittests') if x[0] != '_' and x.endswith('.py')]
