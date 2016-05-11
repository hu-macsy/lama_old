.. _Factory1:

Factory1
========

Factory1 is similiar to :ref:`Factory` but allows one additional argument for the creation of objects.

.. code-block:: c++

  class Base : public Factory1<InputType, OtherType, OutputType>

