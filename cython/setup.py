from setuptools import setup, Extension
from Cython.Build import cythonize
import shutil
import os


sources = ['gintsieve.pyx',
           '../cpp/main.cpp',
           '../cpp/BaseSieve.cpp',
           '../cpp/QuadrantSieve.cpp',
           '../cpp/OctantSieve.cpp',
           '../cpp/DonutSieve.cpp',
           '../cpp/SegmentedSieve.cpp']

extensions = [Extension('gintsieve',
                        sources=sources,
                        extra_compile_args=['-std=c++11', '-stdlib=libc++'],
                        extra_link_args=['-std=c++11', '-stdlib=libc++'],
                        language='c++')]

setup(ext_modules=cythonize(extensions))


# Cleaning up after the cythonize build
shutil.rmtree('build/')

# Rename the shared object file gintsieve.xxxxxx.so
for file in os.listdir('.'):
    if file.endswith('.so'):
        os.rename(file, 'gintsieve.so')
        break
