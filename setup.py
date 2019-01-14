#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

pd2 = Extension(
    'py_power_diagram_2d_double',
    sources = [ 'py_power_diagram/cpp/py_power_diagram_2d_double.cpp' ],
    include_dirs = [ 'ext/pybind11/include' ],
    define_macros = [ ( 'PD_DIM', '2' ), ( 'PD_TYPE', 'double' ), ( 'PD_MODULE_NAME', 'py_power_diagram_2d_double' ) ],
    extra_compile_args = [ '-march=native', '-ffast-math' ],
)

setup(
    name = 'py_power_diagram',
    version = '0.1',
    packages = find_packages( exclude = ( 'tests', 'docs' ) ),
    ext_modules = [ pd2 ],
    install_requires = [
        "numpy",
    ],
)
