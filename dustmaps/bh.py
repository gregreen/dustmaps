#!/usr/bin/env python
#
# bh.py
# Reads the Burstein & Heiles (1982; BH) dust reddening map.
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

import numpy as np
import astropy.coordinates as coordinates
import h5py
import os

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic

def ascii2h5(bh_dir=os.path.join(data_dir_default, 'bh')):
    fname = os.path.join(bh_dir, '{}.ascii')

    f = h5py.File('bh.h5', 'w')

    for region in ('hinorth', 'hisouth'):
        data = np.loadtxt(fname.format(region), dtype='f4')

        # Reshape and clip
        data.shape = (210, 201) # (R, N)
        data = data[:201]   # Last 9 records are empty

        # Use NaNs where no data
        data[data < -9000] = np.nan

        dset = f.create_dataset(
            region,
            data=data,
            chunks=True,
            compression='gzip',
            compression_opts=3
        )

        dset.attrs['axes'] = ('R', 'N')
        dset.attrs['description'] = (
            'HI 21cm column densities, in units of 10*NHYD. '
            'R = 100 + [(90^o-|b|) sin(l)]/[0.3 degrees]. '
            'N = 100 + [(90^o-|b|) cos (l)]/[0.3 degrees].'
        )

    for region in ('rednorth', 'redsouth'):
        data = np.loadtxt(fname.format(region), dtype='f4')

        # Reshape and clip
        data.shape = (94, 1200) # (R, N)
        data = data[:93]   # Last record is empty

        # Use NaNs where no data
        data[data < -9000] = np.nan

        dset = f.create_dataset(
            region,
            data=data,
            chunks=True,
            compression='gzip',
            compression_opts=3
        )

        dset.attrs['axes'] = ('R', 'N')
        dset.attrs['description'] = (
            'E(B-V), in units of 0.001 mag. '
            'R = (|b| - 10) / (0.6 degrees). '
            'N = (l + 0.15) / 0.3 - 1.'
        )

    f.attrs['description'] = (
        'The Burstein & Heiles (1982) dust map.'
    )

    f.close()


class BHQuery(DustMap):
    def __init__(self, bh_dir=os.path.join(data_dir_default, 'bh')):
        f = h5py.File(os.path.join(bh_dir, 'bh.h5'), 'r')
        self._hinorth = f['hinorth'][:]
        self._hisouth = f['hisouth'][:]
        self._rednorth = f['rednorth'][:]
        self._redsouth = f['redsouth'][:]
        f.close()

    def _lb2RN_northcap(self, l, b):
        R = 100. + (90. - b) * np.sin(np.radians(l)) / 0.3
        N = 100. + (90. - b) * np.cos(np.radians(l)) / 0.3
        return np.round(R).astype('i4'), np.round(N).astype('i4')

    def _lb2RN_southcap(self, l, b):
        R = 100. + (90. + b) * np.sin(np.radians(l)) / 0.3
        N = 100. + (90. + b) * np.cos(np.radians(l)) / 0.3
        return np.round(R).astype('i4'), np.round(N).astype('i4')

    def _lb2RN_mid(self, l, b):
        R = (np.abs(b) - 10.) / 0.6
        N = (np.mod(l, 360.) + 0.15) / 0.3 - 1
        return np.round(R).astype('i4'), np.round(N).astype('i4')

    def _lb2ebv_northcap(self, l, b):
        R, N = self._lb2RN_northcap(l, b)
        return -0.0372 + self._hinorth[R,N] * 0.0000357

    def _lb2ebv_southcap(self, l, b):
        R, N = self._lb2RN_southcap(l, b)
        return -0.0372 + self._hisouth[R,N] * 0.0000357

    def _lb2ebv_midnorth(self, l, b):
        R, N = self._lb2RN_mid(l, b)
        return self._rednorth[R,N] * 0.001

    def _lb2ebv_midsouth(self, l, b):
        R, N = self._lb2RN_mid(l, b)
        return self._redsouth[R,N] * 0.001

    @ensure_flat_galactic
    def query(self, coords):
        # gal = coords.transform_to('galactic')
        gal = coords
        l = gal.l.deg
        b = gal.b.deg

        # Detect scalar input
        scalar_input = not hasattr(l, '__len__')
        if scalar_input:
            l = np.array([l])
            b = np.array([b])

        # Fill return array with NaNs
        ebv = np.empty(l.shape, dtype='f8')
        ebv[:] = np.nan

        # Fill northern cap
        idx = (b >= 65.) & (b <= 90.)
        ebv[idx] = self._lb2ebv_northcap(l[idx], b[idx])

        # Fill southern cap
        idx = (b <= -65.) & (b >= -90.)
        ebv[idx] = self._lb2ebv_southcap(l[idx], b[idx])

        # Fill northern midplane
        idx = (b < 65.) & (b >= 10.)
        ebv[idx] = self._lb2ebv_midnorth(l[idx], b[idx])

        # Fill southern midplane
        idx = (b > -65.) & (b <= -10.)
        ebv[idx] = self._lb2ebv_midsouth(l[idx], b[idx])

        if scalar_input:
            ebv = ebv[0]

        return ebv




def main():
    #ascii2h5()
    bh = BHQuery()

    # Calculate E(B-V) on a grid
    l = np.arange(-180, 180., 0.1)
    b = np.arange(-90., 90.01, 0.1)
    l, b = np.meshgrid(l, b)

    c = coordinates.SkyCoord(l, b, frame='galactic', unit='deg')

    ebv = bh.query(c)

    # Apply gamma stretch
    gamma = 0.8
    img = np.power(np.abs(ebv), gamma) * np.sign(ebv)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12,6), dpi=300)
    ax = fig.add_subplot(1,1,1, axisbg='blue')
    ax.imshow(
        img,
        origin='lower',
        interpolation='none',
        cmap='Greys',
        aspect='equal',
        extent=(-180, 180., -90., 90.),
        vmin=-np.power(0.1, gamma),
        vmax=np.power(0.5, gamma),
        rasterized=True
    )
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel(r'$\ell$', fontsize=18)
    ax.set_ylabel(r'$b$', fontsize=18)
    ax.set_title(r'$\mathrm{Burstein - Heiles \ \left( 1982 \right)}$', fontsize=22)

    fig.savefig(os.path.join(output_dir, 'bh.svg'), dpi=300, bbox_inches='tight')
    #plt.show()

    return 0

if __name__ == '__main__':
    main()
