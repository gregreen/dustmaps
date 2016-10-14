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

from __future__ import print_function, division

from setuptools import setup, Extension
from setuptools.command.install import install
import distutils.cmd

import os
import json


class InstallCommand(install):
    description = install.description
    user_options = install.user_options + [
        ('large-data-dir=', None, 'Directory to store large data files in.')
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.large_data_dir = None

    def finalize_options(self):
        if not self.large_data_dir is None:
            self.large_data_dir = os.path.abspath(os.path.expanduser(self.large_data_dir))

        install.finalize_options(self)

    def run(self):
        if not self.large_data_dir is None:
            print('Large data directory is set to: {}'.format(self.large_data_dir))
            with open(os.path.expanduser('~/.dustmapsrc'), 'w') as f:
                json.dump({'data_dir': self.large_data_dir}, f, indent=2)

        # install.do_egg_install(self) # Due to bug in setuptools that causes old-style install
        install.run(self)


def fetch_sfd():
    import dustmaps.sfd
    dustmaps.sfd.fetch()

def fetch_planck():
    import dustmaps.planck
    dustmaps.planck.fetch()

def fetch_bayestar():
    import dustmaps.bayestar
    dustmaps.bayestar.fetch()

def fetch_bh():
    print('Burstein & Heiles (1982) is already installed by default.')


class FetchCommand(distutils.cmd.Command):
    description = ('Fetch dust maps from the web, and store them in the data'
                   'directory.')
    user_options = [
        ('map-name=', None, 'Which map to load.')
    ]

    map_funcs = {
        'sfd': fetch_sfd,
        'planck': fetch_planck,
        'bayestar': fetch_bayestar,
        'bh': fetch_bh
    }

    def initialize_options(self):
        self.map_name = None

    def finalize_options(self):
        try:
            import dustmaps
        except ImportError:
            print('You must install the package dustmaps before running the'
                  'fetch command.')
        if not self.map_name in self.map_funcs:
            print('Valid map names are: {}'.format(self.map_funcs.keys()))

    def run(self):
        print('Fetching map: {}'.format(self.map_name))
        self.map_funcs[self.map_name]()


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='dustmaps',
    version='0.1a2',
    description='Uniform interface for multiple dust reddening maps.',
    long_description=readme(),
    url='https://github.com/gregreen/dustmaps',
    download_url='https://github.com/gregreen/dustmaps/tarball/v0.1a2',
    author='Gregory M. Green',
    author_email='gregorymgreen@gmail.com',
    license='GPLv2',
    packages=['dustmaps'],
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'h5py'
    ],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe=False,
    cmdclass = {
        'install': InstallCommand,
        'fetch': FetchCommand,
    },
)
