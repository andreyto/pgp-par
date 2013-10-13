from distutils.core import setup, Extension

indexSearchExt = Extension('_IndexSearch',
     sources = ['IndexSearch.cpp'])

setup (name = 'IndexSearch',
    version = '1.0',
    description = 'IndexSearch extension build',
    ext_modules = [indexSearchExt])
