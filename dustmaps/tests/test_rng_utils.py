#!/usr/bin/env python
#
# test_rng_utils.py
# Tests for dustmaps.rng_utils.
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

import unittest

import numpy as np

from dustmaps.rng_utils import get_rng


class TestGetRng(unittest.TestCase):
    """The rules for turning a seed_or_rng argument into a generator."""

    def test_none_gives_a_generator(self):
        """No seed gives a usable, unpredictable generator."""
        rng = get_rng()
        self.assertIsInstance(rng, np.random.Generator)

        # Two of them are independent
        other = get_rng()
        self.assertFalse(np.array_equal(rng.integers(0, 10**6, 100),
                                        other.integers(0, 10**6, 100)))

    def test_an_integer_seed_repeats(self):
        """The same integer gives the same draws, however often it is used."""
        first = get_rng(1234).integers(0, 10**6, 100)
        second = get_rng(1234).integers(0, 10**6, 100)
        np.testing.assert_equal(first, second)

    def test_different_integers_differ(self):
        """Different integers give different draws."""
        first = get_rng(1234).integers(0, 10**9, 100)
        second = get_rng(1235).integers(0, 10**9, 100)
        self.assertFalse(np.array_equal(first, second))

    def test_a_generator_is_used_as_it_stands(self):
        """
        A generator is returned unchanged, so that its stream carries on
        instead of restarting.
        """
        rng = np.random.default_rng(7)
        self.assertIs(get_rng(rng), rng)

        first = get_rng(rng).integers(0, 10**9, 100)
        second = get_rng(rng).integers(0, 10**9, 100)
        self.assertFalse(np.array_equal(first, second))

        # ... which is not the same as reseeding it, twice, with its own seed
        third = np.random.default_rng(7).integers(0, 10**9, 100)
        np.testing.assert_equal(first, third)

    def test_a_legacy_random_state_is_refused(self):
        """
        The legacy generator has no ``integers`` method, so it is refused with
        an explanation, rather than raising an obscure error inside a query.
        """
        with self.assertRaises(TypeError) as cm:
            get_rng(np.random.RandomState(7))
        self.assertIn('integers', str(cm.exception))

    def test_a_seed_sequence_is_accepted(self):
        """Anything else is handed to the default generator."""
        seed_seq = np.random.SeedSequence(1234)
        rng = get_rng(seed_seq)
        self.assertIsInstance(rng, np.random.Generator)

        # A SeedSequence, unlike a generator, still gives the same draws
        other = get_rng(np.random.SeedSequence(1234))
        np.testing.assert_equal(rng.integers(0, 10**6, 100),
                                other.integers(0, 10**6, 100))

    def test_a_nonsense_seed_raises(self):
        """A bad seed fails here, rather than producing a stream."""
        with self.assertRaises(TypeError):
            get_rng('birch')


if __name__ == '__main__':
    unittest.main()
