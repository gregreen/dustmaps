dustmaps
========

A uniform interface for a number of 2D and 3D maps of interstellar dust
reddening/extinction.


Supported Dust Maps
-------------------

1. Schlegel, Finkbeiner & Davis (1998; SFD'98)
2. Planck Collaboration (2013)
3. Green, Schlafly, Finbeiner et al. (2015; Bayestar)
4. Burstein & Heiles (1982; BH'82)


Installation
------------

Download the repository and then run:

    python setup.py install --large-data-dir=/path/where/you/want/large/data/files/stored

Alternatively, you can use `pip`:

    pip install dustmaps


Getting the Data
----------------

To fetch the data for the SFD dust map, run:

    python setup.py fetch --map-name=sfd

You can download the other dust maps by changing "sfd" to "planck", "bayestar"
or "bh".

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
    >>> import dustmaps.bayestar
    >>> dustmaps.bayestar.fetch()


Querying the Maps
-----------------

Maps are queried using `astropy.coordinate.SkyCoord` objects. This means that any
coordinate system can be used as input. For example, we can query SFD'98 as
follows:

    >>> from dustmaps.sfd import SFDQuery
    >>> from astropy.coordinates import SkyCoord
    >>>
    >>> sfd = SFDQuery()
    >>>
    >>> c = SkyCoord(
            '05h00m00.00000s',
            '+30d00m00.0000s',
            frame='icrs'
        )
    >>> print sfd(c)
    0.483961

Above, we've used the ICRS coordinate system (the inputs are RA and Dec). We can
use other coordinate systems, such as Galactic coordinates, and we can provide
coordinate arrays. The following example uses both features:

    >>> c = SkyCoord(
            [75.00000000, 130.00000000],
            [-89.00000000, 10.00000000],
            frame='galactic',
            unit='deg'
        )
    >>> print sfd(c)
    [ 0.0146584   0.97695869]
