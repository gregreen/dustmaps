dustmaps
========

A uniform interface for a number of 2D and 3D maps of interstellar dust
reddening/extinction.

Supported Dust Maps
-------------------

1. Schlegel, Finkbeiner & Davis (1998; SFD'98)
2. Burstein & Heiles (1982; BH'82)

Getting the Data
----------------

The data for some dust maps must be downloaded separately.

| Map    | Data included |
| ------ | :-----------: |
| SFD'98 | no            |
| BH'82  | yes           |

The SFD'98 maps can be downloaded [here](http://nebel.rc.fas.harvard.edu/mjuric/lsd-data/sfd-dust-maps/).

Querying the Maps
-----------------

Maps are queried using `astropy.Coordinate` coordinates. This means that any
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
