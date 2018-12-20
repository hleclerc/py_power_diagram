from setuptools import setup, find_packages, Extension

pd2 = Extension(
    'py_power_diagram_2d_double',
    sources = [ 'py_power_diagram/pybind.cpp' ],
    include_dirs = [ 'ext/xsimd/include', 'ext/pybind11/include' ],
    define_macros = [ ( 'PD_DIM', '2' ), ( 'PD_TYPE', 'double' ), ( 'PD_MODULE_NAME', 'py_power_diagram_2d_double' ) ],
    extra_compile_args = [ '-march=native' ]
)
# , '-ffast-math'

setup(
    name = "py_power_diagram",
    version = "0.1",
    packages = [ 'py_power_diagram' ],
    ext_modules = [ pd2 ]
    # packages=find_packages(),
)

# from distutils.core import setup, Extension

# # module1 = Extension('demo',
# #                     define_macros = [('MAJOR_VERSION', '1'),
# #                                      ('MINOR_VERSION', '0')],
# #                     include_dirs = ['/usr/local/include'],
# #                     libraries = ['tcl83'],
# #                     library_dirs = ['/usr/local/lib'],
# #                     sources = ['demo.c'])

# setup(
#     name = 'PowerDiagram',
#     version = '1.0',
#     description = 'This is a demo package',
#     ext_modules = [ module1 ]
# )
