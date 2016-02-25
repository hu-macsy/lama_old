Printable
=========

Instead of the implemation of the operator<< for each class, we use a virtual method within the
class.

.. code-block:: c++

    virtual void writeAt( std::ostream& stream ) const;

All classes derived from the class Printable can override this method to print
info about the actual incarnation.

LAMA has implementated this operator for nearly all classes.

