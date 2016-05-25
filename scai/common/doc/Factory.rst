.. _Factory:

Factory/Factory1
================

A factory provides a static method (here it is always named ``create`` )
that returns an object of a class. But unlike a constructor, the actual object it returns 
might be an instance of a derived class. This factory concept is the central mechanism of the 
LAMA project for any dynamic extension as it allows to create derived objects in software
that is not aware of the extension.

The common library provides the template class ``Factory`` that has two template arguments:

* ``InputType`` is the type used for the choice of the created object
* ``OutputType`` is the type of the created object that is usually a pointer of the class
  for which the factory creates objects.

A (base) class for which a factory is required should derive from this Factory
class. In the following example, a factory becomes available for the class ``Base``, where 
new objects are created by a string value.

.. code-block:: c++

   class Base : public Factory<std::string, Base*>
   {
       ...
   };

Derived classes of Base can register in the factory as follows:

.. code-block:: c++

   class Derived : public Base, Base::Register<Derived>
   {
       static inline std::string createValue()
       {
           return "Derived"; 
       }

       static Base* create()
       {
          return new Derived();
       }
    
       ....
   };

   Base::Register<Derived>::RegisterGuard Base::Register<Derived>::registerGuard;

* By deriving from Base::Register<Derived> registration methods become available
* The variable registerGuard is a guard variable that takes care of the registration 
  in the factory during the static initialization.
* The static method ``createValue`` in the derived class has to define the value (of the
  InputType) used for the creation.
* The static create method will construct the object by using a corresponding constructor
  of the derived class.

For a derived class, objects of this class can now be constructed as follows:

.. code-block:: c++

   Base* obj = Base::create( "Derived" );

In this example, a usual pointer type is used as ``OutputType`` of the factory.
Each call of the routine will create a new object. The caller of the create method takes
the ownership and has to free the object. Therefore it is recommended to use smart pointers.

.. code-block:: c++

   common::unique_ptr<Base> obj ( Base::create( "Derived" ) );

The Factory template class can also be used in such a way that the OutputType is a smart
pointer type. E.g. in the dmemo project (distributed memory) a communicator factory is used
that returns shared pointers for a communicator.

.. code-block:: c++

   class Communicator :: Factory<CommunicatorType, common::shared_ptr<Communicator> > 
   {
      ...
   };

While there is only one factory for the base class, there might be many Registrators, 
usually one for each derived class. Here are some remarks:

* For a derived class that does not derive from the Registrator class in the factory, objects
  of this class cannot be created via the factory.
* Each derived class should use a different create value, i.e. the value returned by the
  static method ``createValue`` should be specific for each derived class.
* It is also possible to use a derived template class.

This is the ways to defined derived template classes that register in the factory for the base class.

.. code-block:: c++

   template<typename T>
   class Derived : public Base, Base::Register<Derived<T> >
   {
   public:

       static inline std::string createValue()
       {
           return typeid( T ).name();
       }
   
       /** Method that creates objects of type Derived2 that will be used for registration. */
   
       static Base* create()
       {
           return new Derived<T>();
       }
    };

Guard for registration should be initiated explicitly.

.. code-block:: c++

    template Base::Register<Derived<float> >::RegisterGuard Base::Register<Derived<float> >::registerGuard;
    template Base::Register<Derived<double> >::RegisterGuard Base::Register<Derived<double> >::registerGuard;

Factory1 is similiar to Factory but allows one additional argument for the creation of objects.

.. code-block:: c++

  class Base : public Factory1<InputType, OtherType, OutputType>

