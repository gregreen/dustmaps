#!/usr/bin/env python
#
# test_override_config.py
# Test code to override default config location
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

import unittest

import os

from ..std_paths import test_dir

class TestConfig(unittest.TestCase):
    def test_config_override(self):
        """
        Test overriding default config location
        """
        os.environ['DUSTMAPS_CONFIG_FNAME'] = os.path.join(test_dir, 'test.json')
        from ..config import config
        self.assertEqual(config['data_dir'], '/my/very/special/path')
        # Remove the entry, but don't overwrite the file
        del config._options['data_dir']
        del os.environ['DUSTMAPS_CONFIG_FNAME']

if __name__ == '__main__':
    unittest.main()
