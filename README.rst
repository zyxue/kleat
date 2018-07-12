======
KLEAT
======


.. image:: https://img.shields.io/pypi/v/kleat.svg
        :target: https://pypi.python.org/pypi/kleat

.. image:: https://img.shields.io/travis/zyxue/kleat.svg
        :target: https://travis-ci.org/zyxue/kleat

.. image:: https://coveralls.io/repos/github/zyxue/kleat/badge.svg
        :target: https://coveralls.io/github/zyxue/kleat

.. image:: https://readthedocs.org/projects/kleat/badge/?version=latest
        :target: https://kleat.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

Cleavage site prediction via de novo assembly.

**Note**: this is a reimplementation from scratch of the following paper (PMID: 25592595_),

.. _25592595: https://www.ncbi.nlm.nih.gov/pubmed/25592595

- Birol I, Raymond A, Chiu R, Nip KM, Jackman SD, Kreitzman M, et al. KLEAT:
  cleavage site analysis of transcriptomes. Pac Symp Biocomput. World
  Scientific. 2015;:347–58.

KLEAT works by

1. Collect polyA evidences (suffix, bridge, link, and blank) and their
   supporting clvs per contig
2. Aggregate polyA evidence per clv identified by (seqname, strand, clv)
3. Find the closest annotated clv per predicted clv and calculate the distance in between.

Suffix was named tail in the original paper, but tail collides with tail of a
bridge contig, which could be confused, so it is renamed to suffix in the code
base.

Blank indicates the contig without any indication of polyA, but could still
support predicting cleavage site since if a transcript is expressed, then a
cleavage site likely exists nearby.

..
   memo: adding hyperlink to a sentence is really awkward in rst!

..
   * Documentation: https://kleat.readthedocs.io.


Install
--------

KLEAT (>3.0.0) supports only Python3 (>=py34). A few key packages include
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

Zero-based index
----------------

Every index is 0-based, including ascii visualization such as

.. code-block::

   Symbols:
   --: ref_skip
   //: hardclip at right
   \\: hardclip at left
   __: deletion
   ┬ : insertion
    └: softclip at left
    ┘: softclip at right

   Abbreviation:
    cc: ctg_clv, clv in contig coordinate
    rc: ref_clv, clv in reference coordinate

   icb: init_clv_beg, initialized beginning index in contig coordinate (for - strand clv)
   irb: init_ref_beg, initialized beginning index in reference coordinate (for - strand clv)

   ice: init_clv_end, initialized end index in contig coordinate (for + strand clv)
   ire: init_ref_end, initialized end index in reference coordinate (for + strand clv)

    TTT
      └AT
    89012 <- one coord (0-based)
      1   <- ten coord

which is different from the display on IGV that is 1-based (although its
underlying system is still 0-based_).

.. _0-based: https://software.broadinstitute.org/software/igv/IGV.


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
