from setuptools import setup, Extension
from Cython.Build import cythonize
import shutil
import os


sources = ['src/gintsieve.pyx',
           '../src/cython_bindings.cpp',
           '../src/BaseSieve.cpp',
           '../src/QuadrantSieve.cpp',
           '../src/OctantSieve.cpp',
           '../src/DonutSieve.cpp',
           '../src/SegmentedSieve.cpp']

extensions = [Extension('gintsieve',
                        sources=sources,
                        extra_compile_args=['-std=c++11', '-stdlib=libc++'],
                        extra_link_args=['-std=c++11', '-stdlib=libc++'],
                        language='c++')]

setup(ext_modules=cythonize(extensions))


# Cleaning up after the cythonize build
shutil.rmtree('build/')
# Removing the cython generated .cpp
os.remove('src/gintsieve.cpp')

# Rename the shared object file gintsieve.xxxxxx.so
for file in os.listdir('.'):
    if file.endswith('.so'):
        os.rename(file, 'gintsieve.so')
        break
