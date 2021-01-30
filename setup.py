"""Install gaussianprimes."""

import os
import setuptools
from Cython.Build import cythonize
import numpy as np


with open('README.md') as f:
  long_description = f.read()
description = 'Generate Gaussian primes.'


sources = [
    'python/gaussianprimes.pyx',
    'src/cython_bindings.cpp',
    'src/BaseSieve.cpp',
    'src/OctantSieve.cpp',
    'src/OctantDonutSieve.cpp',
    'src/SectorSieve.cpp',
    'src/BlockSieve.cpp',
    'src/OctantMoat.cpp'
]

# Calling clang instead of gcc; needed for linux environments
os.environ['CC'] = 'clang'
os.environ['CXX'] = 'clang++'


extensions = [setuptools.Extension(
    'gaussianprimes',
    sources=sources,
    include_dirs=[np.get_include(), 'include'],
    extra_compile_args=['-std=c++11', '-stdlib=libc++'],
    extra_link_args=['-std=c++11', '-stdlib=libc++'],
    language='c++'
)]


setuptools.setup(
    name='gaussianprimes',
    version='1.0.0',
    author='Zeb Engberg',
    author_email='zebengberg@gmail.com',
    url='https://github.com/zebengberg/gaussian-integer-sieve',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',

    python_requires='>=3.7',
    install_requires=['numpy', 'matplotlib', 'setuptools', 'Cython', 'sympy'],
    ext_modules=cythonize(
        extensions,
        compiler_directives={'embedsignature': True}
    ),

    classifiers=['License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Topic :: Scientific/Engineering :: Visualization'],
    license='MIT'
)
