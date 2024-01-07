Install Instructions
====================

What You Need to Get Started
----------------------------

We highly recommend creating a new MJOcast specific conda environment: 

.. code-block:: bash

    conda create --name MJOcastDev python=3.10
    
    
Then install MJOcast via pip: 

MJOcast is available in the [Python Package Index (PyPI)](https://pypi.org/project/MJOcast/). It can be installed using, 
e.g., pip.
    
.. code-block:: bash

    pip3 install MJOcast
    
It can also be installed from the source, which is available on [Zenodo](http://dx.doi.org/10.5281/zenodo.3613752) and [GitHub](https://github.com/willychap/MJOcast). 
Download the source, move into the directory containing the file setup.py and run

.. code-block:: bash

    pip3 install .

The following optional dependency sets can additionaly be installed by adding ["set_name"] behind the above commands:
  * *full_func*: Install these packages to be sure that all options are really available. This might 
    require a higher Python version than for the core functionality alone. RMM may still be computed without these
    additional packages, but the number of alternative approaches is limited.
  * *dev*: packages that are needed for the further development, documentation and testing of MJOcast.

Examples: 
.. code-block:: bash

    pip3 install MJOcast[full_func]

or
.. code-block:: bash

    pip3 install -e .[dev]
    