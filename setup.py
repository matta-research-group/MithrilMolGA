from setuptools import find_packages, setup

# this can be more complicated and clever, but for now it's minimal
setup(
    name='mithrilmolGA',
    version='0.31.0',
    packages=find_packages(include=['mithrilmolGA', 'mithrilmolGA.*']),
)
