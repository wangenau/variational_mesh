.. _installation:

Installation
************

| The code is written for Python 3.4+, but it will work with Python 2.7.
| The following packages are needed for a working installation

* `matplotlib <https://matplotlib.org/>`_
* `numpy <https://numpy.org/>`_
* `pyscf <https://sunqm.github.io/pyscf/>`_

Some examples require extra packages, namely

* `ase <https://wiki.fysik.dtu.dk/ase/>`_
* `pyflosic <https://github.com/pyflosic/pyflosic>`_

The package and all necessary dependencies can be installed using pip

.. code-block:: bash

   $ pip install var-mesh

Alternatively, you can create an installation by downloading the source code

.. code-block:: bash

   $ git clone https://gitlab.com/wangenau/variational_mesh.git
   $ cd variational_mesh
   $ pip install .
