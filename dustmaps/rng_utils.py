#!/usr/bin/env python
#
# rng_utils.py
# Turns a seed or a random number generator into a generator.
#
# Copyright (C) 2026  Gregory M. Green
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

"""
Random number generators for the queries that draw random samples.

The ``query`` functions of the maps that have posterior samples take a
``seed_or_rng`` argument, which the caller may give as an integer, in which case
every call to ``query`` draws the same samples, or as a generator, in which case
the draws carry on from wherever that generator had got to. Sampling modes that
are not random (``mean``, ``median``, ``best``, ...) ignore the argument.

:obj:`get_rng()` does the conversion in one place, so that the query functions
do not each have to know the rules -- in particular that an existing generator
has to be used as it stands, rather than being reseeded (which would defeat the
point of passing one) or handed to :obj:`numpy.random.default_rng` (which does
not accept a generator). The legacy :obj:`numpy.random.RandomState` is refused,
with an explanation, because it has no ``integers`` method, which is what the
query modes draw with.
"""

import numpy as np

__all__ = ['get_rng']


def get_rng(seed_or_rng=None):
    """
    Returns the random number generator to draw from.

    Args:
        seed_or_rng (Optional[:obj:`int`, :obj:`numpy.random.Generator`,
            :obj:`numpy.random.SeedSequence` or
            :obj:`numpy.random.BitGenerator`]): Either a random number
            generator, which is used as it stands, or a seed for a new one.
            Defaults to :obj:`None`, meaning a new generator seeded from the
            operating system, i.e. an unpredictable one.

    Returns:
        (:obj:`numpy.random.Generator`): The generator to draw from. A generator
        that is passed in is returned unchanged, so that its stream continues
        where it left off; anything else is passed to
        :obj:`numpy.random.default_rng`.
    """
    if isinstance(seed_or_rng, np.random.Generator):
        return seed_or_rng

    # The legacy generator is not usable, because the queries draw with
    # `integers`, which it does not have. Say so here, rather than letting an
    # obscure AttributeError surface from inside a query.
    if isinstance(seed_or_rng, np.random.RandomState):
        raise TypeError(
            'seed_or_rng cannot be a numpy.random.RandomState, because it has '
            'no `integers` method. Pass an integer, or a generator such as '
            'numpy.random.default_rng(seed).')

    return np.random.default_rng(seed_or_rng)
