.. _troubleshooting:

Trouble Shooting
----------------

Here we list known problems when building LAMA and their solutions.

ABI Error
^^^^^^^^^

When you use C11 capable GCC and you get error messages containing GLIBCXX and ABI in it you can use

.. code-block:: bash

	_GLIBCXX_USE_CXX11_ABI

as described in this post on <a href="http://stackoverflow.com/questions/34571583/understanding-gcc-5s-glibcxx-use-cxx11-abi-or-the-new-abi" target="_blank">stackoverflow</a>.
