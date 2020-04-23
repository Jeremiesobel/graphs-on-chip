#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="graphs-on-chip",
    version="0.1-dev",
    description="Analyze confocal image stacks",
    author="Jeremie Sobel, Gustave Ronteix",
    author_email="gronteix@pasteur.fr",
    url="https://github.com/microfluidix/graphs-on-chip",
    install_requires=[
        "numpy",
        "scipy",
        "networkx",
        "pytest",
        "tqdm",
        "pandas",
        "click"
    ],
    python_requires=">=3.8",
    packages=find_packages(),
)
