.. role:: python(code)
   :language: python

Examples
========

Getting Started
---------------

Here, we'll look up the reddening at a number of different locations on the sky.
We specify coordinates on the sky using
`astropy.coordinates.SkyCoord <http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_
objects. This allows us a great deal of flexibility in how we specify sky
coordinates. We can use different coordinate frames (e.g.,
`Galactic <https://en.wikipedia.org/wiki/Galactic_coordinate_system>`_,
`equatorial <https://en.wikipedia.org/wiki/Equatorial_coordinate_system>`_,
`ecliptic <https://en.wikipedia.org/wiki/Ecliptic_coordinate_system>`_),
different units (e.g., degrees, radians,
`hour angles <https://en.wikipedia.org/wiki/Hour_angle>`_), and either
scalar or vector input.

For our first example, let's load the
`Schlegel, Finkbeiner & Davis (1998) <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_
-- or "SFD" -- dust reddening map, and then query the reddening at one location
on the sky:

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.sfd import SFDQuery
    
    coords = SkyCoord('12h30m25.3s', '15d15m58.1s', frame='icrs')
    sfd = SFDQuery()
    ebv = sfd(coords)

    print('E(B-V) = {:.3f} mag'.format(ebv))
    
    >>> E(B-V) = 0.030 mag

A couple of things to note here:

1. In this example, we used :python:`from __future__ import print_function` in
   order to ensure compatibility with both Python 2 and 3.
2. Above, we used the
   `ICRS coordinate system <https://en.wikipedia.org/wiki/International_Celestial_Reference_System>`_,
   by specifying :python:`frame='icrs'`.
3. :python:`SFDQuery` returns reddening in a unit that is similar to magnitudes
   of **E(B-V)**. However, care should be taken: a unit of SFD reddening is not
   quite equivalent to a magnitude of **E(B-V)**. The way to correctly convert
   SFD units to extinction in various broadband filters is to use the
   conversions in
   `Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_.

We can query the other maps in the :code:`dustmaps` package with only minor
modification to the above code. For example, here's how we would query the
Planck Collaboration (2013) dust map:

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    
    coords = SkyCoord('12h30m25.3s', '15d15m58.1s', frame='icrs')
    planck = PlanckQuery()
    ebv = planck(coords)
    
    print('E(B-V) = {:.3f} mag'.format(ebv))
    
    >>> E(B-V) = 0.035 mag


Querying Reddening at an Array of Coordinates
---------------------------------------------

We can also query an array of coordinates, as follows:


.. code-block :: python
    
    from __future__ import print_function
    import numpy as np
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    from dustmaps.sfd import SFDQuery
    
    l = np.array([0., 90., 180.])
    b = np.array([15., 0., -15.])
    
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    
    planck = PlanckQuery()
    planck(coords)
    >>> array([ 0.50170666,  1.62469053,  0.29259142])
    
    sfd = SFDQuery()
    sfd(coords)
    >>> array([ 0.55669367,  2.60569382,  0.37351534], dtype=float32)

The input need not be a flat array. It can have any shape -- the shape of the
output will match the shape of the input:

.. code-block :: python
    
    from __future__ import print_function
    import numpy as np
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    
    l = np.linspace(0., 180., 12)
    b = np.zeros(12, dtype='f8')
    l.shape = (3, 4)
    b.shape = (3, 4)
    
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    
    planck = PlanckQuery()
    
    ebv = planck(coords)
    
    print(ebv)
    >>> [[ 315.52438354   28.11778831   23.53047562   20.72829247]
         [   2.20861101   15.68559361    1.46233201    1.70338535]
         [   0.94013882    1.11140835    0.38023439    0.81017196]]
    
    print(ebv.shape)
    >>> (3, 4)


Querying 3D Reddening Maps
--------------------------

When querying a 3D dust map, there are two slight complications:

1. There is an extra axis -- distance -- to care about.
2. Many 3D dust maps are probabilistic, so we need to specify whether we want
   the median reddening, mean reddening, a random sample of the reddening, etc.

Let's see how this works out with the "Bayestar" dust map of
`Green, Schlafly & Finkbeiner (2015) <http://argonaut.skymaps.info>`_.

How Distances are Handled
~~~~~~~~~~~~~~~~~~~~~~~~~

If we don't provide distances in our input, :code:`dustmaps` will assume we want dust
reddening along the entire line of sight.

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.bayestar import BayestarQuery
    
    coords = SkyCoord(180., 0., unit='deg', frame='galactic')
    
    # Note that below, we could use version='bayestar2017' to get the newer
    # version of the map. Note, however, that the reddening units are not
    # identical in the two versions of the map. See Green et al. (2018) for
    # an explanation of the units.
    bayestar = BayestarQuery(max_samples=2, version='bayestar2015')
    
    ebv = bayestar(coords, mode='random_sample')
    
    print(ebv)
    >>> [ 0.00476     0.00616     0.0073      0.00773     0.00796     0.07453
          0.07473     0.0748      0.07807     0.07831     0.18957999  0.2013
          0.20448001  0.20734     0.21008     0.73733997  0.75415999  0.93702
          0.93956     1.09001005  1.09141004  1.11407995  1.11925006  1.12212002
          1.12284994  1.12289     1.12296999  1.12305999  1.12308002  1.12309003
          1.12311995]

Here, the Bayestar map has given us a single random sample of the cumulative
dust reddening *along the entire line of sight* -- that is, to a set of
distances. To see what those distances are, we can call:

.. code-block :: python
    
    bayestar.distances
    >>> <Quantity [  0.06309573,  0.07943282,  0.1       ,  0.12589255,
                     0.15848933,  0.19952621,  0.25118864,  0.31622776,
                     0.3981072 ,  0.50118726,  0.63095725,  0.79432821,
                     1.        ,  1.2589252 ,  1.58489335,  1.99526215,
                     2.51188707,  3.1622777 ,  3.98107076,  5.01187277,
                     6.3095727 ,  7.94328403, 10.        , 12.58925152,
                    15.84893322, 19.95262146, 25.11886978, 31.62277603,
                    39.81070709, 50.11872864, 63.09572601] kpc>

The return type is an `astropy.unit.Quantity <http://astropy.readthedocs.io/en/stable/api/astropy.units.Quantity.html>`_
instance, which keeps track of units.

If we provide Bayestar with distances, then it will do the distance
interpolation for us, returning the cumulative dust reddening out to specific
distances:

.. code-block :: python
    
    import astropy.units as units
    
    coords = SkyCoord(180.*units.deg, 0.*units.deg,
                      distance=500.*units.pc, frame='galactic')
    ebv = bayestar(coords, mode='median')
    
    print(ebv)
    >>> 0.10705789

Because we have explicitly told Bayestar what distance to evaluate the map at,
it returns only a single value.


How Probability is Handled
~~~~~~~~~~~~~~~~~~~~~~~~~~

The Bayestar 3D dust map is probabilistic, meaning that it stores random samples
of how dust reddening could increase along each sightline. Sometimes we might be
interested in the median reddening to a given point in space, or we might want
to have all the samples of reddening out to that point. We specify how we want
to deal with the probabilistic nature of the map by providing the keyword
argument :code:`mode` to :code:`dustmaps.bayestar.BayestarQuery.__call__`.

For example, if we want all the reddening samples, we invoke:

.. code-block :: python
    
    l = np.array([30.,  60., 90.]) * units.deg
    b = np.array([10., -10., 15.]) * units.deg
    d = np.array([1.5,  0.3, 4.0]) * units.kpc
    
    coords = SkyCoord(l, b, distance=d, frame='galactic')
    
    ebv = bayestar(coords, mode='samples')
    
    print(ebv.shape) # (# of coordinates, # of samples)
    >>> (3, 2)
    
    print(ebv)
    >>> [[ 0.24641787  0.27142054]    # Two samples at the first coordinate
         [ 0.01696703  0.0149225 ]    # Two samples at the second coordinate
         [ 0.08348     0.11068   ]]   # Two samples at the third coordinate

If we instead ask for the mean reddening, the shape of the output is different:

.. code-block :: python
    
    ebv = bayestar(coords, mode='mean')
    
    print(ebv.shape) # (# of coordinates)
    >>> (3,)
    
    print(ebv)
    >>> [ 0.25891921  0.09121627  0.09708   ]

The only axis is for the different coordinates, because we have reduced the
samples axis by taking the mean.

In general, the shape of the output from the Bayestar map is:

.. code-block :: python
    
    (coordinate, distance, sample)

where any of the axes can be missing (e.g., if only one coordinate was
specified, if distances were provided, or if the median reddening was
requested).

Percentiles are handled in much the same way as samples. In the following
query, we request the 16th, 50th and 84th percentiles of reddening at each
coordinate, using the same coordinates as we generated in the previous example:

.. code-block :: python
    
    ebv = bayestar(coords, mode='percentile', pct=[16., 50., 84.])
    
    print(ebv)
    >>> [[ 0.24789949  0.25583497  0.26986977]  # Percentiles at 1st coordinate
         [ 0.01505572  0.01814967  0.02750403]  # Percentiles at 2nd coordinate
         [ 0.0860716   0.09787634  0.10787529]] # Percentiles at 3rd coordinate

We can also pass a single percentile:

.. code-block :: python
    
    ebv = bayestar(coords, mode='percentile', pct=25.)
    
    print(ebv)
    >>> [ 0.24930404  0.01524667  0.08961   ] # 25th percentile at 3 coordinates


Getting Quality Assurance Flags from the Bayestar Dust Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the Bayestar dust maps, one can retrieve QA flags by providing the keyword
argument :code:`return_flags=True`:

.. code-block :: python
    
    ebv, flags = bayestar(coords, mode='median', return_flags=True)
    
    print(flags.dtype)
    >>> [('converged', '?'), ('reliable_dist', '?')]
    
    print(flags['converged']) # Whether or not fit converged in each pixel
    >>> [ True  True  True]
    
    # Whether or not map is reliable at requested distances
    print(flags['reliable_dist'])
    >>> [ True False  True]

If the coordinates do not include distances, then instead of
:code:`'reliable_dist'`, the query will return the minimum and maxmimum reliable
distance moduli of the map in each requested coordinate:

.. code-block :: python
    
    l = np.array([30.,  60., 90.]) * units.deg
    b = np.array([10., -10., 15.]) * units.deg
    
    coords = SkyCoord(l, b, frame='galactic')
    
    ebv, flags = bayestar(coords, mode='median', return_flags=True)
    
    print(flags['min_reliable_distmod'])
    >>> [ 7.875       8.24800014  6.87300014]
    
    print(flags['max_reliable_distmod'])
    >>> [ 15.18599987  15.25500011  15.00699997]

We can see from the above that in the previous example, the reason the second
coordinate was labeled unreliable was because the requested distance (300 pc)
was closer than a distance modulus of 8.248 (corresponding to ~450 pc).


Plotting the Dust Maps
----------------------

We'll finish by plotting a comparison of the SFD, Planck Collaboration and
Bayestar Dust maps. First, we'll import the necessary modules:

.. code-block :: python
    
    from __future__ import print_function
    
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    
    import astropy.units as units
    from astropy.coordinates import SkyCoord
    
    from dustmaps.sfd import SFDQuery
    from dustmaps.planck import PlanckQuery
    from dustmaps.bayestar import BayestarQuery

Next, we'll set up a grid of coordinates to plot, centered on the Aquila South
cloud:

.. code-block :: python
    
    l0, b0 = (37., -16.)
    l = np.arange(l0 - 5., l0 + 5., 0.05)
    b = np.arange(b0 - 5., b0 + 5., 0.05)
    l, b = np.meshgrid(l, b)
    coords = SkyCoord(l*units.deg, b*units.deg,
                      distance=1.*units.kpc, frame='galactic')

Then, we'll load up and query three different dust maps:

.. code-block :: python
    
    sfd = SFDQuery()
    Av_sfd = 2.742 * sfd(coords)
    
    planck = PlanckQuery()
    Av_planck = 3.1 * planck(coords)
    
    bayestar = BayestarQuery(max_samples=1)
    Av_bayestar = 2.742 * bayestar(coords)

We've assumed :math:`R_V = 3.1`, and used the coefficient from
`Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_
to convert SFD and Bayestar reddenings to magnitudes of :math:`A_V`.

Finally, we create the figure using :code:`matplotlib`:

.. code-block :: python
    
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

Here's the result:

.. image :: figs/comparison.png


Querying the web server
-----------------------

Some of the maps included in this package are large, and can take up a lot
of memory, or be slow to load. To make it easier to work with these maps,
some of them are available to query over the internet. As of now, the
following maps can be queried remotely:

* Bayestar (all versions)
* SFD

The API for querying these maps remotely is almost identical to the API
for local queries. For example, the following code queries SFD remotely:

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.sfd import SFDWebQuery
    
    l = [180., 160.]
    b = [30., 45.]
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    sfd = SFDWebQuery()
    ebv = sfd(coords)
    
    print(ebv)
    
    >>> [0.04704102 0.02022794]

The following example queries the Bayestar2019 dust map remotely. The web
interface takes the same arguments as the local interface:

.. code-block :: python
    
    import astropy.units as u
    from dustmaps.bayestar import BayestarWebQuery
    
    l = [90., 150., 35.] * u.deg
    b = [10., 12., -25.] * u.deg
    d = [500., 3500., 1000.] * u.pc
    coords = SkyCoord(l, b, distance=d, frame='galactic')
    
    q = BayestarWebQuery(version='bayestar2019')
    E = q(coords, mode='median')
    
    print(E)
    
    >>> [0.13       0.63       0.09999999]

The :code:`query_gal()` and :code:`query_equ()` convenience functions also
work with web queries. Continuing from the previous example,

.. code-block :: python
    
    E = q.query_gal([120., 125.], [-5., -10.],
                    d=[1.5, 1.3],
                    mode='random_sample')
    print(E)
    
    >>> [0.32 0.24]

Please take it easy on our web server. If you want to query multiple
coordinates, then bundle them up into one query. If you want to query
a *very large* number of coordinates, consider downloading the maps and
querying them locally instead.
