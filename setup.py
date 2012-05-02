# -*- encoding: utf8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("eigsym3d", ["eigsym3d.pyx","src/eigsym3d.c"], extra_compile_args=['-ffast-math']),
]

setup(
    name = "pyEigSym3D",
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    author = 'Alexis Mignon',
	author_email = "alexis.mignon@gmail.com",
    url = 'http://github.com/alexis-mignon/pyEigSym3D',
    description = 'Fast eigen decomposition of 3x3 symmetric matrices',
    license = 'GPL 2.0',
)

