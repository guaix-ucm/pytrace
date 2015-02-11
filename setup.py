
from __future__ import print_function

from setuptools import setup, Extension
from setuptools import find_packages

import sys

setup(name='pytrace',
      version='0.1.dev',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      license='GPLv3',
      description='Tracing experiment',
      packages=find_packages('.'),
      long_description=open('README.md').read()
      )
