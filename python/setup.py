from setuptools import setup, Extension
from Cython.Build import cythonize
import shutil
import os
import numpy as np


sources = ['src/gintsieve.pyx',
           '../src/cython_bindings.cpp',
           '../src/BaseSieve.cpp',
           '../src/QuadrantSieve.cpp',
           '../src/OctantSieve.cpp',
           '../src/DonutSieve.cpp',
           '../src/SegmentedSieve.cpp']

extensions = [Extension('gintsieve',
                        sources=sources,
                        include_dirs=[np.get_include()],
                        extra_compile_args=['-std=c++11', '-stdlib=libc++'],
                        extra_link_args=['-std=c++11', '-stdlib=libc++'],
                        language='c++')]

setup(
    name='gintsieve',
    version='1.0',
    url='https://github.com/zebengberg/gaussian-integer-sieve',
    author='Zeb Engberg',
    license='MIT',
    ext_modules=cythonize(extensions))


try:
    print('Cleaning up after build...')
    # Cleaning up build directory after the cythonize build
    shutil.rmtree('build/')
    # Removing the cython generated .cpp
    os.remove('src/gintsieve.cpp')
except FileNotFoundError:
    pass

# Rename the shared object file gintsieve.xxxxxx.so
for file in os.listdir('.'):
    if file.endswith('.so'):
        os.rename(file, 'gintsieve.so')
        break
