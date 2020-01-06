from setuptools import setup, Extension
from Cython.Build import cythonize

import os
# instead of clang try gcc
# os.environ["CC"] = "g++-9"



extensions = [Extension('gintsieve',
                        sources=['gintsieve.pyx'],
                        extra_compile_args=['-std=c++11', '-stdlib=libc++'],
                        extra_link_args=['-std=c++11', '-stdlib=libc++'],
                        language='c++'
                        )]


setup(ext_modules=cythonize(extensions))


