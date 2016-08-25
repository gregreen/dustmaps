#!/usr/bin/env python
#
# setup.py
# Package "dustmaps" for pip.
#
# Copyright (C) 2016  Gregory M. Green
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from setuptools import setup, Extension

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='dustmaps',
    version='0.1a0',
    description='Uniform interface for multiple dust reddening maps.',
    long_description=readme(),
    url='https://github.com/gregreen/dustmaps',
    author='Gregory M. Green',
    author_email='gregorymgreen@gmail.com',
    license='GPLv2',
    packages=['dustmaps'],
    install_requires=[
        'numpy',
        'scipy',
        'astropy'
        'h5py'
    ],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe=False
)
