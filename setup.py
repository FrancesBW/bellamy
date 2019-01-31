#! /usr/bin/env python
"""
Setup for MACCAS
"""
import os
import sys
from setuptools import setup

reqs=['numpy>=1.15.4', 'astropy>=2.0, <3','scipy>=1.1.0','argparse>=1.2.1', 'psutil>=5.4.8', 'matplotlib>=2.2.3']

setup(
    name="maccas",
    version="0.1",
    author="Frances Buckland-Willis",
    description="MACCAS: Multi Attribute Cross-matcher with Correction for wArped Sky",
    url="https://github.com/FrancesBW/radio_cross_matching",
    long_description=['README.md'],
    packages=['functions'],
    install_requires=reqs,
    scripts=['scripts/maccas'],
    classifiers=['Programming Language :: Python :: 2.7']
)
