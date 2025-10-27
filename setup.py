from setuptools import find_packages, setup

# Minimal setup matching the folder name
setup(
    name='MithrilMolGA',  # Use the exact folder name (case-sensitive)
    version='0.1.0',
    packages=find_packages(include=['MithrilMolGA', 'MithrilMolGA.*']),
)

