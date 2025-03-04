#!/usr/bin/env python
#
# plot_bayestar.py
# An example of how to query the "Bayestar" dust map of
# Green, Schlafly, Finkbeiner et al. (2015).
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

from dustmaps.bayestar import BayestarQuery


def numpy2pil(a, vmin, vmax):
    a = np.clip((a - vmin) / (vmax - vmin), 0., 1.)
    a = (254.99 * a).astype('u1')
    return PIL.Image.fromarray(a)


def main():
    w,h = (2056,1024)
    l_0 = 130.

    # Set up Bayestar query object
    print('Loading bayestar map...')
    bayestar = BayestarQuery(max_samples=1)

    # Create a grid of coordinates
    print('Creating grid of coordinates...')
    l = np.linspace(-180.+l_0, 180.+l_0, 2*w)
    b = np.linspace(-90., 90., 2*h+2)
    b = b[1:-1]
    l,b = np.meshgrid(l, b)

    l += (np.random.random(l.shape) - 0.5) * 360./(2.*w)
    b += (np.random.random(l.shape) - 0.5) * 180./(2.*h)

    ebv = np.empty(l.shape+(3,), dtype='f8')

    for k,d in enumerate([0.5, 1.5, 5.]):
        # d = 5.    # We'll query integrated reddening to a distance of 5 kpc
        coords = SkyCoord(l*u.deg, b*u.deg, d*u.kpc, frame='galactic')

        # Get the dust median reddening at each coordinate
        print('Querying map...')
        ebv[:,:,k] = bayestar.query(coords, mode='median')

    ebv[:,:,2] -= ebv[:,:,1]
    ebv[:,:,1] -= ebv[:,:,0]

    # Convert the output array to a PIL image and save
    print('Saving image...')
    img = numpy2pil(ebv[::-1,::-1,:], 0., 1.5)
    img = img.resize((w,h), resample=PIL.Image.LANCZOS)
    fname = 'bayestar.png'
    img.save(fname)

    return 0


if __name__ == '__main__':
    main()
