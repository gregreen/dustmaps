#!/usr/bin/env python
#
# dustexceptions.py
# Defines exceptions for the dustmaps package.
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

from . import std_paths

class Error(Exception):
    pass

class CoordFrameError(Error):
    pass


def data_missing_message(package, name):
    return ("The {name} dust map is not in the data directory:\n\n"
            "    {data_dir}\n\n"
            "To change the data directory, call:\n\n"
            "    from dustmaps.config import config\n"
            "    config['data_dir'] = '/path/to/data/directory'\n\n"
            "To download the {name} map to the data directory, call:\n\n"
            "    import dustmaps.{package}\n"
            "    dustmaps.{package}.fetch()\n").format(
                data_dir=std_paths.data_dir(),
                package=package,
                name=name)
