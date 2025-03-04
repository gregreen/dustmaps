#!/usr/bin/env python
#
# plot_chen.py
# An example of how to query the Chen et al. (2014) 3D dust map.
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

from __future__ import print_function

import numpy as np
import os.path

try:
    import PIL.Image
except ImportError as error:
    print('This example requires Pillow or PIL.\n'
          'See <http://pillow.readthedocs.io/en/stable/installation.html>.')
    raise error

from astropy.coordinates import SkyCoord
import astropy.units as u

from dustmaps.chen2014 import Chen2014Query


def numpy2pil(a, vmin, vmax, fill=0):
    mask = np.isnan(a)
    a = np.clip((a - vmin) / (vmax - vmin), 0., 1.)
    a = (254.99 * a).astype('u1')
    a[mask] = fill
    return PIL.Image.fromarray(a)


def main():
    w,h = (2056/2, 2056/2)

    # Set up Chen2014 object
    print('Loading Chen+(2014) map...')
    query = Chen2014Query()

    # Create a grid of coordinates
    print('Creating grid of coordinates...')
    # l = np.linspace(140., 240., 2*w)
    # b = np.linspace(-60., 40., 2*h)
    l = np.linspace(186., 202., 2*w)
    b = np.linspace(-24., -8., 2*h)
    dl = l[1] - l[0]
    db = b[1] - b[0]
    l,b = np.meshgrid(l, b)

    # l += (np.random.random(l.shape) - 0.5) * dl
    # b += (np.random.random(l.shape) - 0.5) * db

    A = np.empty(l.shape+(3,), dtype='f8')

    for k,d in enumerate([0.5, 1.0, 4.]):
        coords = SkyCoord(l*u.deg, b*u.deg, d*u.kpc, frame='galactic')

        # Get the dust median reddening at each coordinate
        print('Querying map to {:.1f} kpc...'.format(d))
        A[:,:,k] = query.query(coords)

        # Convert the output array to a PIL image and save
        print('Saving image...')
        img = numpy2pil(A[::-1,::-1,k], 0., 2., fill=255)
        img = img.resize((w,h), resample=PIL.Image.LANCZOS)
        fname = 'chen2014_{:03.1f}kpc.png'.format(d)
        img.save(fname)

    A[:,:,2] -= A[:,:,1]
    A[:,:,1] -= A[:,:,0]

    # Convert the output array to a PIL image and save
    print('Saving image...')
    img = numpy2pil(A[::-1,::-1,:], 0., 2., fill=255)
    img = img.resize((w,h), resample=PIL.Image.LANCZOS)
    fname = 'chen2014.png'
    img.save(fname)

    return 0


if __name__ == '__main__':
    main()
