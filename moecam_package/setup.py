#!/usr/bin/env python3
"""
MOECAM Package Setup
===================

Multi-Objective Evolutionary Comparison and Analysis Module
"""

from setuptools import setup, find_packages

# Package setup
setup(
    name="moecam",
    version="1.0.0",
    description="Multi-Objective Evolutionary Comparison and Analysis Module",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "cffi>=1.14.0",
    ],
    setup_requires=[
        "cffi>=1.14.0",
        "pybind11>=2.6"
    ],
    cffi_modules=["moecam/cffi_build.py:ffibuilder"],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
)
