#!/usr/bin/env python
#
# plot_iphas.py
# An example of how to query the "IPHAS" dust map of
# Sale et al. (2014).
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

from dustmaps.iphas import IPHASQuery


def numpy2pil(a, vmin, vmax, fill=0):
    mask = np.isnan(a)
    a = np.clip((a - vmin) / (vmax - vmin), 0., 1.)
    a = (254.99 * a).astype('u1')
    a[mask] = fill
    return PIL.Image.fromarray(a)


def main():
    w,h = (2*2056, 2*int(2056*(30./200.)))
    l_0 = 122.5

    # Set up IPHASquery object
    print('Loading IPHAS map...')
    iphas = IPHASQuery()

    # Create a grid of coordinates
    print('Creating grid of coordinates...')
    l = np.linspace(-100.+l_0, 100.+l_0, 2*w)
    b = np.linspace(-15., 15., 2*h)
    dl = l[1] - l[0]
    db = b[1] - b[0]
    l,b = np.meshgrid(l, b)

    l += (np.random.random(l.shape) - 0.5) * dl
    b += (np.random.random(l.shape) - 0.5) * db

    A = np.empty(l.shape+(3,), dtype='f8')

    for k,d in enumerate([0.5, 1.5, 5.]):
        # d = 5.    # We'll query integrated reddening to a distance of 5 kpc
        coords = SkyCoord(l*u.deg, b*u.deg, d*u.kpc, frame='galactic')

        # Get the dust median reddening at each coordinate
        print('Querying map...')
        A[:,:,k] = iphas.query(coords, mode='random_sample')

    A[:,:,2] -= A[:,:,1]
    A[:,:,1] -= A[:,:,0]

    # Convert the output array to a PIL image and save
    print('Saving image...')
    img = numpy2pil(A[::-1,::-1,:], 0., 4.5, fill=255)
    img = img.resize((w,h), resample=PIL.Image.LANCZOS)
    fname = 'iphas.png'
    img.save(fname)

    return 0


if __name__ == '__main__':
    main()
