Installation
============

There are two ways to install :code:`dustmaps`.


1. Using :code:`pip`
--------------------

From the commandline, run

.. code-block :: bash

    pip install dustmaps

You may have to use :code:`sudo`.

Next, we'll configure the package and download the dust maps we'll want to use.
Start up a python interpreter and type:

.. code-block :: python

    from dustmaps.config import config
    config['data_dir'] = '/path/to/store/maps/in'

    import dustmaps.sfd
    dustmaps.sfd.fetch()

    import dustmaps.planck
    dustmaps.planck.fetch()

    import dustmaps.bayestar
    dustmaps.bayestar.fetch()

    import dustmaps.iphas
    dustmaps.iphas.fetch()

    import dustmaps.marshall
    dustmaps.marshall.fetch()

    import dustmaps.chen2014
    dustmaps.chen2014.fetch()

All the dust maps should now be in the path you gave to
:code:`config['data_dir']`. Note that these dust maps can be very large - some
are several Gigabytes! Only download those you think you'll need.


2. Using :code:`setup.py`
-------------------------

An alternative way to download :code:`dustmaps`, if you don't want to use
:code:`pip`, is to download or clone the respository from
https://github.com/gregreen/dustmaps. Then, from the root directory of the
package, run

.. code-block :: bash

    python setup.py install --large-data-dir=/path/to/store/maps/in

Then, fetch the maps you'd like to use:

.. code-block :: bash

    python setup.py fetch --map-name=sfd
    python setup.py fetch --map-name=planck
    python setup.py fetch --map-name=bayestar
    python setup.py fetch --map-name=iphas
    python setup.py fetch --map-name=marshall
    python setup.py fetch --map-name=chen2014

Since these maps are very large - up to several Gigabytes - be careful to only
download those you think you'll need. That's it!
