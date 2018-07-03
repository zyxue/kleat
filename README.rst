======
kleat
======


.. image:: https://img.shields.io/pypi/v/kleat.svg
        :target: https://pypi.python.org/pypi/kleat

.. image:: https://img.shields.io/travis/zyxue/kleat.svg
        :target: https://travis-ci.org/zyxue/kleat

.. image:: https://readthedocs.org/projects/kleat/badge/?version=latest
        :target: https://kleat.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Cleavage site prediction via de novo assembly

..
   * Documentation: https://kleat.readthedocs.io.


Install
--------

kleat (>3.0.0) supports only Python3 (>=py34). A few key packages include
pysam_, pandas_, scikit-learn_.

.. _pysam: https://github.com/pysam-developers/pysam
.. _pandas: https://github.com/pandas-dev/pandas
.. _scikit-learn: https://github.com/scikit-learn/scikit-learn

First, it's recommended to create a virtual environment, using either
conda_

.. _conda: https://conda.io/miniconda.html

.. code-block:: bash

   conda create -p venv-kleat python=3
   source activate venv-kleat/

or pip_ + virtualenv_

.. _pip: https://github.com/pypa/pip
.. _virtualenv: https://github.com/pypa/virtualenv

.. code-block:: bash

   pip install virtualenv # skip this step if virtualenv is available already
   virtualenv venv-kleat
   . venv-kleat/bin/activate

Then install kleat with pip

.. code-block:: bash

   pip install git+https://github.com/zyxue/kleat.git#egg=kleat


Features
--------

* TODO


Development
-----------

.. code-block:: bash

   virtualenv venv
   . venv/bin/activate
   pip install -r requirements_dev.txt
   python setup.py develop
   kleat --help

To uninstall

.. code-block:: bash

   python setup.py develop --uninstall


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
