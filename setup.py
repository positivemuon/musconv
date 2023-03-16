#!/usr/bin/env python
##https://gist.github.com/juancarlospaco/75d258d6ffb6e15f395c
from setuptools import find_packages, setup
import json 


if __name__ == '__main__':
    with open('setup.json','r') as info:
        kwargs = json.load(info)
    setup(
        packages=find_packages(),
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        **kwargs
    )

