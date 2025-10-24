#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup script for pymhm package
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "Python package for mesoscale Hydrological Model"

# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="pymhm",
    version="0.1.0",
    author="Sanjeev Bashyal",
    author_email="sanjeev.bashyal01@gmail.com",
    description="Python package for mesoscale Hydrological Model",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/SanjeevBashyal/pymhm",
    project_urls={
        "Bug Tracker": "https://github.com/SanjeevBashyal/pymhm/issues",
        "Documentation": "https://github.com/SanjeevBashyal/pymhm",
        "Source Code": "https://github.com/SanjeevBashyal/pymhm",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Hydrology",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    include_package_data=True,
    package_data={
        "pymhm": [
            "*.ui",
            "*.qrc",
            "*.png",
            "*.txt",
            "i18n/*.ts",
            "i18n/*.qm",
        ],
    },
    entry_points={
        "console_scripts": [
            "pymhm=pymhm.cli:main",
        ],
    },
    keywords="hydrology, modeling, mesoscale, water, environment, qgis, plugin",
    zip_safe=False,
)
