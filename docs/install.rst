Installation
============

GFFtk can be installed using pip, conda, or from source.

Using pip
---------

.. code-block:: bash

    pip install gfftk

Using conda
-----------

GFFtk is available through the bioconda channel:

.. code-block:: bash

    conda install -c bioconda gfftk

From source
-----------

To install from source, clone the repository and install using pip:

.. code-block:: bash

    git clone https://github.com/nextgenusfs/gfftk.git
    cd gfftk
    pip install -e .

Dependencies
-----------

GFFtk requires the following Python packages:

* natsort
* numpy
* requests
* gb-io

These dependencies will be automatically installed when using pip or conda.

Verifying Installation
---------------------

To verify that GFFtk has been installed correctly, run:

.. code-block:: bash

    gfftk --version

This should display the version number of GFFtk.

Updating GFFtk
-------------

To update GFFtk to the latest version:

Using pip:

.. code-block:: bash

    pip install --upgrade gfftk

Using conda:

.. code-block:: bash

    conda update -c bioconda gfftk
