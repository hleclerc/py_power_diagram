#!/usr/bin/env python3

from setuptools import setup, find_packages

with open( 'requirements.txt' ) as f:
    requirements = f.read().splitlines()
    
setup(
    name = 'py_power_diagram',
    version = '0.1',
    install_requires = requirements,
    packages = [ 'py_power_diagram' ]
)
