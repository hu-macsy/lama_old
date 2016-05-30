Class ContextFunction
=====================

.. code-block:: c++

  template<typename FunctionType>
  class ContextFunction : public _ContextFunction
  {
  public:
      ContextFunction(); /** Only default Constructor at this time */
      ContextFunction( const ContextFunction& other );
      FunctionType get( ContextType ctx ) const; 
      void set( ContextType ctx, FunctionType fn );
      using _ContextFunction::validContext;
      ContextType validContext( ContextType preferedCtx );
      ContextType validContext( const _ContextFunction& other, ContextType preferedCtx );
  };

The class ``ContextFunction`` has the following properties:

- an object of the class ``ContextFunction`` is just an array of typed function pointers

- the array of typed function pointers is indexed by the context type (enum value)

- base class ``_ContextFunction`` is used as an array of untyped function pointers (for internal purpose only)

- validContext returns 

  - preferred context if function pointer has been set 

  - otherwise the first context type for which it has been registered (usually Host).