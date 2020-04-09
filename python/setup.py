from setuptools import setup, Extension
from Cython.Build import cythonize
import shutil
import os
import numpy as np


sources = ['src/gintsieve.pyx',
           '../src/cython_bindings.cpp',
           '../src/BaseSieve.cpp',
           '../src/OctantSieve.cpp',
           '../src/OctantDonutSieve.cpp',
           '../src/SectorSieve.cpp',
           '../src/BlockSieve.cpp',
           '../src/OctantMoat.cpp']

# Calling clang instead of gcc; needed for linux environments
os.environ['CC'] = 'clang'
os.environ['CXX'] = 'clang++'

extensions = [Extension('gintsieve',
                        sources=sources,
                        include_dirs=[np.get_include()],
                        # last flag below disables depreciated np warning
                        extra_compile_args=['-std=c++11', '-stdlib=libc++', '-Wno-#warnings'],
                        extra_link_args=['-std=c++11', '-stdlib=libc++'],
                        language='c++')]

print('\nInstalling...')
setup(
    name='gintsieve',
    version='1.0',
    url='https://github.com/zebengberg/gaussian-integer-sieve',
    author='Zeb Engberg',
    author_email='zebengberg@gmail.com',
    license='MIT',
    install_requires = ['numpy', 'matplotlib', 'setuptools'],
    ext_modules=cythonize(extensions))



print('\nCleaning up after build...')
# Rename the shared object file from gintsieve.xxxxxx.so to gintsieve.so
for file in os.listdir('.'):
    if file.endswith('.so'):
        os.rename(file, 'gintsieve.so')
        break

# Cleaning up build directory after the cythonize build
shutil.rmtree('build/')
# Removing the cython generated .cpp
os.remove('src/gintsieve.cpp')


print('\nTesting the module just built...')
import gintsieve
try:
    # Values taken from http://oeis.org/A091100
    assert gintsieve.count_gprimes(10 ** 0) == 0
    assert gintsieve.count_gprimes(10 ** 1) == 16
    assert gintsieve.count_gprimes(10 ** 2) == 100
    assert gintsieve.count_gprimes(10 ** 3) == 668
    assert gintsieve.count_gprimes(10 ** 4) == 4928
    assert gintsieve.count_gprimes(10 ** 5) == 38404
    assert gintsieve.count_gprimes(10 ** 6) == 313752
    assert gintsieve.count_gprimes(10 ** 7) == 2658344
    assert gintsieve.count_gprimes(10 ** 8) == 23046512
    assert gintsieve.count_gprimes(10 ** 9) == 203394764
except AssertionError:
    print('Tests failed. Something went very wrong! The module gintsieve is corrupt.\n')
else:
    print('Tests passed. Module gintsieve is ready to go.\n')
