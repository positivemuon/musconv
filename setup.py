# coding: utf-8
##https://gist.github.com/juancarlospaco/75d258d6ffb6e15f395c
from setuptools import find_packages, setup
from json import load as json_load

if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        setup_kwargs = json_load(info)
    setup(
        packages=find_packages(exclude=['docs', 'tests','examples']),
        **setup_kwargs
    )

