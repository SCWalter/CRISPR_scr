#!/usr/bin/env python

from setuptools import setup, find_packages
setup(
    name = "crisprscr",
    version = "0.1",
    packages = find_packages(),
    scripts = ['say_hello.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'hello': ['*.msg'],
    },
    
    entry_points={
        'console_scripts': [
            'crispr-scr = my_package.some_module:main_func',
        ]
    }

    # metadata for upload to PyPI
    author = "Yunhai Luo",
    author_email = "me@example.com",
    description = "This is a pipeline for analyzing CRISPR screen related data",
    license = "MIT",
    keywords = "CRISPR screen pipeline",
    url = "https://github.com/SCWalter/CRISPR_scr",   # project home page

    # could also include long_description, download_url, classifiers, etc.
)
)
