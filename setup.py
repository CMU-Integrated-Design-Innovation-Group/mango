from setuptools import setup

with open('requirements.txt') as f:
    required_packages = f.read().splitlines()

setup(
    install_requires=required_packages
)