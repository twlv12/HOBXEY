from setuptools import setup
from Cython.Build import cythonize

setup(
    name='hobxey',
    ext_modules=cythonize("main.pyx"),
    zip_safe=False,
)