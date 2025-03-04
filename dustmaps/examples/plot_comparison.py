
#!/usr/bin/env python
#
# plot_comparison.py
# An example of how to plot three different dust maps.
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as units
from astropy.coordinates import SkyCoord

from dustmaps.sfd import SFDQuery
from dustmaps.planck import PlanckQuery
from dustmaps.bayestar import BayestarQuery

import os


def main():
    c0 = SkyCoord.from_name('orion a', frame='galactic')
    print(c0)

    # l = np.arange(c0.l.deg - 5., c0.l.deg + 5., 0.05)
    # b = np.arange(c0.b.deg - 5., c0.b.deg + 5., 0.05)

    l0, b0 = (37., -16.)
    l = np.arange(l0 - 5., l0 + 5., 0.05)
    b = np.arange(b0 - 5., b0 + 5., 0.05)
    l, b = np.meshgrid(l, b)
    coords = SkyCoord(l*units.deg, b*units.deg,
                      distance=1.*units.kpc, frame='galactic')

    sfd = SFDQuery()
    Av_sfd = 2.742 * sfd(coords)

    planck = PlanckQuery()
    Av_planck = 3.1 * planck(coords)

    bayestar = BayestarQuery(max_samples=1)
    Av_bayestar = 2.742 * bayestar(coords)

    fig = plt.figure(figsize=(12,4), dpi=150)

    for k,(Av,title) in enumerate([(Av_sfd, 'SFD'),
                                   (Av_planck, 'Planck'),
                                   (Av_bayestar, 'Bayestar')]):
        ax = fig.add_subplot(1,3,k+1)
        ax.imshow(
            np.sqrt(Av)[::,::-1],
            vmin=0.,
            vmax=2.,
            origin='lower',
            interpolation='nearest',
            cmap='binary',
            aspect='equal'
        )
        ax.axis('off')
        ax.set_title(title)

    fig.subplots_adjust(wspace=0., hspace=0.)
    plt.savefig('comparison.png', dpi=150)

    return 0


if __name__ == '__main__':
    main()
