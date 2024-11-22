#!/usr/bin/env python
#
# std_paths.py
# Defines a set of paths used by scripts in the dustmaps module.
#
# Copyright (C) 2016  Gregory M. Green
#
# dustmaps is free software: you can redistribute it and/or modify
# it under the terms of either:
#
# - The GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version, or
# - The 2-Clause BSD License (also known as the Simplified BSD License).
#
# You should have received copies of the GNU General Public License
# and the BSD License along with this program.
#

from __future__ import print_function, division

import os
from .config import config


script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir_default = os.path.abspath(os.path.join(script_dir, 'data'))
test_dir = os.path.abspath(os.path.join(script_dir, 'tests'))
output_dir_default = os.path.abspath(os.path.join(script_dir, 'output'))


def fix_path(path):
    """
    Returns an absolute path, expanding both '~' (to the user's home
    directory) and other environmental variables in the path.
    """
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))


def data_dir():
    """
    Returns the directory used to store large data files (e.g., dust maps).
    """
    dirname = config.get('data_dir', data_dir_default)
    return fix_path(dirname)


def output_dir():
    """
    Returns a directory that can be used to store temporary output.
    """
    dirname = config.get('output_dir', output_dir_default)
    return fix_path(dirname)
