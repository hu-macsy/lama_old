.. _TypeList:

TypeList
========

In C++, a type list can be defined by a recursive template pattern like the following one:

.. code-block:: c++

    template<typename H, typename T>
    struct TypeList
    {
        typedef H head;
        typedef T tail;
    };

A NullType is used for termination.

.. code-block:: c++

   class NullType {};

Here is an example:

.. code-block:: c++

   TypeList<int, TypeList<float, TypeList<double, NullType> > >

.. code-block:: c++

   #define SCAI_TYPELIST( int, float, dobule )

In the following we show how a type list can be used to call a template method for a certain 
list of types.

.. code-block:: c++

   template<typename T>
   void output()
   {
       ...
   }

General defintion of Calling for template list:

.. code-block:: c++

   template<typename TList> struct Calling;

Here is  the termination call

.. code-block:: c++

   template<> struct Calling<mepr::NullType>
   {
       static void call()
       {
       }
   };

Call output for header T and recursive call for tail of list

.. code-block:: c++

   template<typename HeadType, typename TailTypes>
   struct Calling<mepr::TypeList<HeadType, TailTypes> >
   {
       static void call()
       {
           output<HeadType>();
           Calling<TailTypes>::call();
       }
   };


