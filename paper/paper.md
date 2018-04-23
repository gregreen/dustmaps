---
title: 'dustmaps: A Python interface for maps of interstellar dust'
tags:
  - Python
  - astronomy
  - interstellar medium
  - interstellar reddening
  - interstellar extinction
authors:
  - name: Gregory M. Green
    orcid: 0000-0001-5417-2260
    affiliation: 1
affiliations:
  - name: Porat Fellow, Stanford/KIPAC
    index: 1
date: 19 April 2018
bibliography: paper.bib
---

# Summary

Correcting for interstellar dust extinction is a critical step in many analyses of astrophysical data. Indeed, a standard dust reddening map, @Schlegel:1998, is one of the highest cited papers in astrophysics.

The ``dustmaps`` package provides a uniform Python interface for several commonly used maps of interstellar dust, including two-dimensional maps such as @Schlegel:1998, @Planck:2013 and @Lenz:2017, and three-dimensional maps such as @Marshall:2006 and @Green:2015. ``dustmaps`` makes use of ``Astropy``'s coordinate-system package [``astropy.coordinates.SkyCoord``, @astropy], making it easy to query dust maps in a wide variety of coordinate systems (Equatorial, Galactic, Ecliptic, etc.). Additionally, ``dustmaps`` handles the downloading of the supported dust maps for users, and allows users to query some dust maps from a remote server, avoiding the need to download large data files.

Development of ``dustmaps`` takes place on GitHub [@github_dustmaps], and any issues with the software or feature suggestions (e.g., the addition of new dust maps) should be raised there.

An example of the type of analysis which can be carried out with ``dustmaps`` is given below. The left panel is a plot of dust reddening in @Green:2018 to a distance of 800 pc, while the right panel shows the correlation between @Green:2018 and @Planck:2013.

![Example of the type of analysis made easy by ``dustmaps``.](figure.pdf)

# References
