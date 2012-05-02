# -*- encoding: utf8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("eigsym3d", ["eigsym3d.pyx","src/eigsym3d.c"], extra_compile_args=['-ffast-math']),
]

setup(
    name = "eigsym3d",
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

