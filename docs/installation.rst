Installation
============

There are two ways to install :code:`dustmaps`.

1. Using :code:`pip`
--------------------

From the commandline, run

.. code-block :: bash

    pip install dustmaps

You may have to use :code:`sudo`.

Next, we'll configure the package and download the dust maps we'll want to use. Start up a python interpreter and type:::

    from dustmaps.config import config
    config['data_dir'] = '/path/to/store/maps/in'

    import dustmaps.sfd
    dustmaps.sfd.fetch()

    import dustmaps.planck
    dustmaps.planck.fetch()

    import dustmaps.bayestar
    dustmaps.bayestar.fetch()

All the dust maps should now be in the path you gave to :code:`config['data_dir']`. Note that these dust maps can be very large - some are several Gigabytes! Only download those you think you'll need.

2. Using :code:`setup.py`
-------------------------

First, download the package from https://github.com/gregreen/dustmaps. Then, from the root directory of the package, run

.. code-block :: bash

    python setup.py install --large-data-dir=/path/to/store/maps/in

Then, fetch the maps you'd like to use:

.. code-block :: bash

    python setup.py fetch --map-name=sfd
    python setup.py fetch --map-name=planck
    python setup.py fetch --map-name=bayestar

Since these maps are very large - up to several Gigabytes - be careful to only download those you think you'll need. That's it!
