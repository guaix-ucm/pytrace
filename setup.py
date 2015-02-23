
from __future__ import print_function

from setuptools import setup, Extension
from setuptools import find_packages

import sys

# try to handle gracefully Cython
try:
    from Cython.Distutils import build_ext
    ext1 = Extension('trace._traces',
                     ['trace/traces.pyx', 'trace/fitter.cpp',
                      'trace/Trace.cpp'],
                     language='c++')
    cmdclass = {'build_ext': build_ext}
except ImportError:
    print('We do not have Cython, just using the generated files')
    ext1 = Extension('trace._traces',
                     ['trace/traces.cpp', 'trace/fitter.cpp',
                      'trace/Trace.cpp'],
                     language='c++')
    cmdclass = {}

setup(name='pytrace',
      version='0.1.dev',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      license='GPLv3',
      description='Tracing experiment',
      packages=find_packages('.'),
      ext_modules=[ext1],
      cmdclass=cmdclass,
      long_description=open('README.md').read()
      )
