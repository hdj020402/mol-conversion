#!/usr/bin/env python3
"""Setup script for mol-conversion package"""

from setuptools import setup, find_packages

# Simple configuration without external dependencies
setup(
    name="mol_conversion",
    version="0.1.0",
    description="Molecular conversion tools for XYZ format processing",
    author="Dejun Hu",
    author_email="hudejun2002@gmail.com",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/hdj020402/mol-conversion",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "rdkit>=2020.03.5",
        "tqdm>=4.60.0",
        # Note: openbabel must be installed via conda: conda install openbabel -c conda-forge
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=22.0",
            "flake8>=4.0",
        ],
        # openbabel must be installed via conda: conda install openbabel -c conda-forge
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    include_package_data=True,
)