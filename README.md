[![DOI](http://joss.theoj.org/papers/10.21105/joss.00695/status.svg)](https://doi.org/10.21105/joss.00695) [![DOI](https://zenodo.org/badge/59614814.svg)](https://zenodo.org/badge/latestdoi/59614814)

dustmaps
========

The ``dustmaps`` package provides a uniform interface for dealing with a number
of 2D and 3D maps of interstellar dust reddening/extinction.

Supported Dust Maps
-------------------

The currently supported dust maps are:

1. Burstein & Heiles (1982; BH'82)
2. Chen et al. (2014)
3. Green, Schlafly, Finkbeiner et al. (2015,2018,2019; Bayestar)
4. Marshall et al. (2006)
5. Planck Collaboration (2013)
6. Planck Collaboration (2016; GNILC)
7. Sale et al. (2014; IPHAS)
8. Schlegel, Finkbeiner & Davis (1998; SFD'98)
9. Lenz, Hensley & Doré (2017)
10. Peek & Graves (2010)
11. Leike & Enßlin (2019)
12. Leike, Glatzle & Enßlin (2020)

To request addition of another dust map in this package, [file an issue on
GitHub](https://github.com/gregreen/dustmaps/issues), or submit a pull request.


Installation
------------

Download the repository from [GitHub](https://github.com/gregreen/dustmaps) and
then run:

    python setup.py install --large-data-dir=/path/where/you/want/large/data/files/stored

Alternatively, you can use the Python package manager `pip`:

    pip install dustmaps


Getting the Data
----------------

To fetch the data for the SFD dust map, run:

    python setup.py fetch --map-name=sfd

You can download the other dust maps by changing "sfd" to "planck",
"planckGNILC", "bayestar", "iphas", "marshall", "chen2014", "lenz2017",
"pg2010", "leikeensslin2019", "leike2020" or "bh".

Alternatively, if you have used `pip` to install `dustmaps`, then you can
configure the data directory and download the data by opening up a python
interpreter and running:

    >>> from dustmaps.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>>
    >>> import dustmaps.sfd
    >>> dustmaps.sfd.fetch()
    >>>
    >>> import dustmaps.planck
    >>> dustmaps.planck.fetch()
    >>>
    >>> import dustmaps.planck
    >>> dustmaps.planck.fetch(which='GNILC')
    >>>
    >>> import dustmaps.bayestar
    >>> dustmaps.bayestar.fetch()
    >>>
    >>> import dustmaps.iphas
    >>> dustmaps.iphas.fetch()
    >>>
    >>> import dustmaps.marshall
    >>> dustmaps.marshall.fetch()
    >>>
    >>> import dustmaps.chen2014
    >>> dustmaps.chen2014.fetch()
    >>>
    >>> import dustmaps.lenz2017
    >>> dustmaps.lenz2017.fetch()
    >>>
    >>> import dustmaps.pg2010
    >>> dustmaps.pg2010.fetch()
    >>>
    >>> import dustmaps.leike_ensslin_2019
    >>> dustmaps.leike_ensslin_2019.fetch()
    >>>
    >>> import dustmaps.leike2020
    >>> dustmaps.leike2020.fetch()


Querying the Maps
-----------------

Maps are queried using
[`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord)
objects. This means that any coordinate system supported by `astropy` can be
used as input. For example, we can query SFD'98 as follows:

    >>> from dustmaps.sfd import SFDQuery
    >>> from astropy.coordinates import SkyCoord
    >>>
    >>> sfd = SFDQuery()
    >>>
    >>> c = SkyCoord(
            '05h00m00.00000s',
            '+30d00m00.0000s',
            frame='icrs')
    >>> print sfd(c)
    0.483961

Above, we have used the ICRS coordinate system (the inputs are RA and Dec). We
can use other coordinate systems, such as Galactic coordinates, and we can
provide coordinate arrays. The following example uses both features:

    >>> c = SkyCoord(
            [75.00000000, 130.00000000],
            [-89.00000000, 10.00000000],
            frame='galactic',
            unit='deg')
    >>> print sfd(c)
    [ 0.0146584   0.97695869]


Documentation
-------------

Read the full documentation at http://dustmaps.readthedocs.io/en/latest/.


Citation
--------

If you make use of this software in a publication, please cite
[Green (2018) in The Journal of Open Source Software](https://doi.org/10.21105/joss.00695):

    @ARTICLE{2018JOSS....3..695M,
           author = {{Green}, {Gregory M.}},
            title = "{dustmaps: A Python interface for maps of interstellar dust}",
          journal = {The Journal of Open Source Software},
             year = "2018",
            month = "Jun",
           volume = {3},
           number = {26},
            pages = {695},
              doi = {10.21105/joss.00695},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2018JOSS....3..695M},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }


Development
-----------

Development of `dustmaps` takes place on GitHub, at
https://github.com/gregreen/dustmaps. Any bugs, feature requests, pull requests,
or other issues can be filed there. Contributions to the software are welcome.
