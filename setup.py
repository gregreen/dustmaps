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
import io


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

def fetch_bayestar(**kwargs):
    import dustmaps.bayestar
    dustmaps.bayestar.fetch(**kwargs)

def fetch_iphas():
    import dustmaps.iphas
    dustmaps.iphas.fetch()

def fetch_marshall():
    import dustmaps.marshall
    dustmaps.marshall.fetch()

def fetch_chen2014():
    import dustmaps.chen2014
    dustmaps.chen2014.fetch()

def fetch_leikeensslin2019():
    import dustmaps.leike_ensslin_2019
    dustmaps.leike_ensslin_2019.fetch()

def fetch_leike2020():
    import dustmaps.leike2020
    dustmaps.leike2020.fetch()

def fetch_lenz2017():
    import dustmaps.lenz2017
    dustmaps.lenz2017.fetch()

def fetch_pg2010():
    import dustmaps.pg2010
    dustmaps.pg2010.fetch()

def fetch_bh():
    print('Burstein & Heiles (1982) is already installed by default.')


class FetchCommand(distutils.cmd.Command):
    description = ('Fetch dust maps from the web, and store them in the data '
                   'directory.')
    user_options = [
        ('map-name=', None, 'Which map to load.')]

    map_funcs = {
        'sfd': fetch_sfd,
        'planck': fetch_planck,
        'bayestar': fetch_bayestar,
        'bayestar2015': lambda: fetch_bayestar(version='bayestar2015'),
        'bayestar2017': lambda: fetch_bayestar(version='bayestar2017'),
        'bayestar2019': lambda: fetch_bayestar(version='bayestar2019'),
        'bh': fetch_bh,
        'iphas': fetch_iphas,
        'marshall': fetch_marshall,
        'chen2014': fetch_chen2014,
        'lenz2017': fetch_lenz2017,
        'pg2010': fetch_pg2010,
        'leikeensslin2019': fetch_leikeensslin2019,
        'leike2020': fetch_leike2020
    }

    def initialize_options(self):
        self.map_name = None

    def finalize_options(self):
        try:
            import dustmaps
        except ImportError:
            print('You must install the package dustmaps before running the '
                  'fetch command.')
        if not self.map_name in self.map_funcs:
            print('Valid map names are: {}'.format(self.map_funcs.keys()))

    def run(self):
        print('Fetching map: {}'.format(self.map_name))
        self.map_funcs[self.map_name]()


def readme():
    with io.open('README.md', mode='r', encoding='utf-8') as f:
        return f.read()


setup(
    name='dustmaps',
    version='1.0.6',
    description='Uniform interface for multiple dust reddening maps.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/gregreen/dustmaps',
    download_url='https://github.com/gregreen/dustmaps/archive/v1.0.6.tar.gz',
    author='Gregory M. Green',
    author_email='gregorymgreen@gmail.com',
    license='GPLv2',
    packages=['dustmaps'],
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'h5py',
        'healpy',
        'requests',
        'progressbar2',
        'six'
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
