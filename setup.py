from setuptools import find_packages, setup

setup(
    name='astromol',
    packages=find_packages(),
    package_data={'astromol' : ['bibtex/*.bib']},
    version='2021.0.0',
    description='',
    author='Brett A. McGuire (Department of Chemistry, Massachusetts Institute of Technology)',
    license='MIT',
)