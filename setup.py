#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" setup"""
import json

from setuptools import find_packages, setup

if __name__ == "__main__":
    with open("setup.json", "r", encoding="utf-8") as info:
        kwargs = json.load(info)

    with open("README.md", encoding="utf-8") as file:
        l_descr = file.read()

    setup(
        packages=find_packages(),
        long_description=l_descr,
        long_description_content_type="text/markdown",
        **kwargs
    )
