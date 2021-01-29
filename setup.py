from setuptools import setup, Extension

import shutil
import os
import numpy as np


# try:
#   print('lala')
from Cython.Build import cythonize
# except ModuleNotFoundError:

#   # create closure for deferred import

#   def cythonize(*args, ** kwargs):
#     print('ok')
#     from Cython.Build import cythonize
#     return cythonize(*args, ** kwargs)


sources = [
    'python/gintsieve.pyx',
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


extensions = [Extension(
    'gintsieve',
    sources=sources,
    include_dirs=[np.get_include(), 'include'],
    # last flag below disables depreciated np warning
    extra_compile_args=['-std=c++11', '-stdlib=libc++'],  # '-Wno-#warnings'],
    extra_link_args=['-std=c++11', '-stdlib=libc++'],
    language='c++'
)]

print('\nInstalling...')
setup(
    name='gintsieve',
    version='1.0.0',
    url='https://github.com/zebengberg/gaussian-integer-sieve',
    author='Zeb Engberg',
    author_email='zebengberg@gmail.com',
    license='MIT',
    python_requires='>=3',
    install_requires=['numpy', 'matplotlib', 'setuptools', 'Cython'],
    ext_modules=cythonize(extensions)
    # ext_modules=extensions
)


# print('\nCleaning up after build...')
# # Rename the shared object file from gintsieve.xxxxxx.so to gintsieve.so
# for file in os.listdir('.'):
#   if file.endswith('.so'):
#     os.rename(file, 'gintsieve.so')
#     break

# # Cleaning up build directory after the cythonize build
# shutil.rmtree('build/')
# # Removing the cython generated .cpp
# os.remove('src/gintsieve.cpp')
