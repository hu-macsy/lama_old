Enum Types
==========

* should always be defined in an own namespace (unless they are defined in struct or class)
* operator<< should always be defined in same namespace

The operator<< for enum types might be used in logging messages and should print useful names
instead of simple int values.

.. code-block:: c++

    namespace context
    {
        enum ContextType
        {
            Host, 
            CUDA,
            OpenCL,
            MIC, 
            UserContext,
            MaxContext
        };

        std::ostream& operator<<( std::ostream& stream, const ContextType& type );

    } 

Due to the own namespace name clashes are avoided, each enum value must be qualified via the namespace.

.. code-block:: c++

    context::ContextType ctype = context::CUDA;

