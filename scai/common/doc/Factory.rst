Factory 
-------

LAMA makes like many other object-oriented software heavy use of inheritance. Typically for each
base class there are many derived classes available.

Wanted: applications that can deal with all kind of derived classes even with dynamically added ones (by library module)

.. code-block:: c++

   ContextPtr cudaContext = Context::create( "CUDA" );
   shared_ptr<_Harray> myArray ( _HArray::create( ScalarType::FLOAT ) );

A Factory is a class that allows the dynamic creation of derived objects for a base class without
knowing at compile time which derived classes might be available.

.. code-block:: c++

   class Base : public Factory<InputType, OutputType>


.. code-block:: c++

   class Base : public Factory<std::string, Base*>, public Printable
   {
      public:
      ...
   };


.. code-block:: c++

  vector<string> values;  // string is create type for the factory

  Base::getCreateValues( values );

  cout << "Factory of Base: " << endl;

  for ( size_t i = 0; i < values.size(); ++i )
  {
      cout << "   Registered values[" << i << "] = " << values[i] << endl;
  }

  Base* obj1 = Base::create( "D1" );
  Base* obj2 = Base::create( "D2" );

  if ( Base::canCreate( "D3" ) )
  {
      ...
  }

This is the ways to defined derived classes that register in the factory for the base class.

.. code-block:: c++

    class Derived1 : public Base, Base::Register<Derived1>
    {
    public:

        /** This static functions provides the value for which this class registers in the factory. */

        static inline std::string createValue()
        {
            return "D1";
        }

        /** Method that creates objects of type Derived2 that will be used for registration. */
    
        static inline Base* create()
        {
            return new Derived1();
        }
    
        ...
    };

A guard for the registration should be instantiated explicitly in the code.

.. code-block:: c++

    Base::Register<Derived1>::RegisterGuard Base::Register<Derived1>::register

For template classes this mechanism is as follows:

.. code-block:: c++

   template<typename T>
   class Derived : public Base, Base::Register<Derived<T> >
   {
   public:

       /** Provide the value for which this class registers in the factory. */

       static inline std::string createValue()
       {
           return typeid( T ).name();
       }

       /** Method that creates objects, will be used for registration. */

       static Base* create()
       {
           return new Derived<T>();
       }
       ...
   };

Guard for registration should be initiated explicitly.

.. code-block:: c++

   Base::Register<Derived<float> >::RegisterGuard Base::Register<Derived<float> >::register;
   Base::Register<Derived<double> >::RegisterGuard Base::Register<Derived<double> >::register;

Factory1 is similiar to Factory but allows one additional argument for the creation of objects.

.. code-block:: c++

  class Base : public Factory1<InputType, OtherType, OutputType>

