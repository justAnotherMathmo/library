from setuptools import setup

__version__ = '1.0.0'

with open("README", 'r') as f:
    long_description = f.read()

setup(
    name='library',
    version=__version__,
    description='Common functions, mainly for Project Euler problems',
    long_description=long_description,
    author='justAnotherMathmo',
    author_email='kchap21@hotmail.com',
    packages=['library'],
    install_requires=['numpy'],
)
