from distutils.core import setup, Extension
from Cython.Build import cythonize

import os
# instead of clang try gcc
# os.environ["CC"] = "g++-9"



extensions = [Extension('gintsieve',
                        sources=['gintsieve.pyx'],
                        extra_compile_args=['-std=c++11', '-stdlib=libc++', '-v']
                        )]


setup(ext_modules=cythonize(extensions))


