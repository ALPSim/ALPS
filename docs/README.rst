Dependencies
============
Compilation of the documentation depends on the Sphinx, Breathe and doxygen. You can install Sphinx and Breathe on any machine (including Windows) by using the Pythin Index Package manager ``pip``

.. code-block:: bash
   
   pip install sphinx
   pip install breathe

You can also use your systems native package manager in many cases.

Compilation
===========
In order to compile the newest documentation run:

.. code-block:: bash

   make html

The generated documentation can be found in ``_build/html/``. You can also generate other formats such as PDF.
