
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

#include <iostream>
#include <unistd.h>

using namespace scai;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "MICExample" )

template<typename ValueType>
ValueType sum( const ValueType array[], const IndexType n )
{
    ValueType zero = static_cast<ValueType>( 0 );

    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )

    ValueType val = 0;

    const void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n ), out( val )
    {
        val = 0;

        const ValueType* array = static_cast<const ValueType*>( arrayPtr );

        #pragma omp parallel for reduction( +:val )

        for ( IndexType i = 0; i < n; ++i )
        {
            val += array[i];
        }
    }

    return val;
}

template<typename ValueType>
void add( ValueType* array, const IndexType n )
{
    SCAI_LOG_INFO( logger, "add # array = " << array << ", n = " << n )

    void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );

        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] += 1;
        }
    }
}

void printContextFactory()
{
    std::vector<context::ContextType> values;

    Context::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::cout << "Registered values[" << i << "] = " << values[i] << std::endl;
    }
}

int main()
{
    printContextFactory();

    std::cout << "try to get " << context::MIC << " context from factory" << std::endl;
    ContextPtr micContext = Context::getContextPtr( context::MIC );
    SCAI_ASSERT( micContext, "NULL context" )
    std::cout << "micContext = " << *micContext << std::endl;

    MemoryPtr micMemory = micContext->getMemoryPtr();
    SCAI_ASSERT( micMemory, "NULL memory" )
    std::cout << "micMemory = " << *micMemory << std::endl;

    std::cout << "try to get " << context::Host << " context from factory" << std::endl;
    ContextPtr hostContext = Context::getContextPtr( context::Host );
    SCAI_ASSERT( hostContext, "NULL context" )
    std::cout << "hostContext = " << *hostContext << std::endl;

    MemoryPtr hostMemory = hostContext->getMemoryPtr();
    SCAI_ASSERT( hostMemory, "NULL memory" )
    std::cout << "hostMemory = " << *hostMemory << std::endl;

    MemoryPtr micHostMemory = micContext->getHostMemoryPtr();
    SCAI_ASSERT( micHostMemory, "NULL memory" )
    std::cout << "micHostMemory = " << *micHostMemory << std::endl;

    const IndexType N = 5;

    LAMAArray<double> data( micContext );
    
    std::cout << "data = " << data << std::endl;

    {
        SCAI_LOG_INFO( logger, "write only on host" )

        WriteOnlyAccess<double> writeData( data, N );

        double* dataHost = writeData;

        for ( IndexType i = 0; i < N; ++i )
        {
            writeData[i] = 1.0;
        }
    }

    std::cout << "After host write: data = " << data << std::endl;

    {
        SCAI_LOG_INFO( logger, "read on mic" )
        ReadAccess<double> read( data, micContext );
        SCAI_CONTEXT_ACCESS( micContext )
        double s = sum( read.get(), data.size() );
        std::cout << "sum = " << s << ", should be " << N  << std::endl;
    }

    std::cout << "After mic read: data = " << data << std::endl;

    {
        SCAI_LOG_INFO( logger, "write on mic" )
        WriteAccess<double> write( data, micContext );
        SCAI_CONTEXT_ACCESS( micContext )
        add( static_cast<double*>( write ), data.size() );
    }

    std::cout << "After mic write: data = " << data << std::endl;

    {
        SCAI_LOG_INFO( logger, "read on host" )
        ReadAccess<double> read( data );
        sleep( 1 );
        for ( IndexType i = 0; i < N; ++i )
        {
            SCAI_ASSERT_EQUAL( read[i], 2 * 1.0, "wrong value after add, i = " << i )
        }
    }

    std::cout << "After host read: data = " << data << std::endl;
}

