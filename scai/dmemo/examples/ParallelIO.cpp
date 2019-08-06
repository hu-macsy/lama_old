
#include <scai/hmemo.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CollectiveFile.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>

using namespace scai;
using namespace dmemo;

#define HOST_PRINT( rank, msg )                 \
    {                                           \
        if ( rank == 0 )                        \
        {                                       \
            std::cout << msg << std::endl;      \
        }                                       \
    }                                           \

template<typename ValueType>
static hmemo::HArray<ValueType> distributedArray( const dmemo::Distribution& dist, ValueType ( *fill )( IndexType ) )
{
    hmemo::HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {
        IndexType localIndex = 0;   // running local index

        for ( auto& entry : hmemo::hostWriteOnlyAccess( localArray, dist.getLocalSize() ) )
        {
            entry = fill( dist.local2Global( localIndex++ ) );
        }

    }  // filled the local array with 'global' values

    return localArray;    // each processor gets its local part
}

typedef int ValueType;

int main( int argc, const char* argv[] )
{
    if ( argc < 2 )
    {
        std::cout << "Usage: " << argv[0] << " filename                       // read the file" << std::endl;
        std::cout << "Usage: " << argv[0] << " filename  <nentries> [repeat]  // write the file" << std::endl;
        return -1;
    }

    const char* fileName = argv[1];

    bool out = ( argc > 2 );

    auto fillArray = []( IndexType k ) { return ValueType( 2 * k + 1 ); };

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    auto rank = comm->getRank();

    if ( out )
    { 
        SCAI_REGION( "main.out" )

        const IndexType N = atoi( argv[2] );
        const IndexType K = argc > 3 ? atoi( argv[3] ) : 1;

        HOST_PRINT( rank, "Write collective file " << fileName << " with " << K << " x " << N << " entries" )

        auto dist = blockDistribution( N, comm );
    
        hmemo::HArray<ValueType> myData = distributedArray<ValueType>( *dist, fillArray );
    
        auto time = common::Walltime::get();

        auto outFile = comm->collectiveFile();
        outFile->open( fileName, "w" );
        outFile->writeSingle( K );
 
        for ( IndexType iter = 0; iter < K; iter ++ )
        {
            SCAI_REGION( "main.write" )
            outFile->writeSingle( N );
            outFile->writeAll( myData, dist->lb() );
        }

        auto nBytes = outFile->getOffset();

        outFile->close();

        time = common::Walltime::get() - time;
        
        double rateGBs = static_cast<double>( nBytes ) / ( 1000.0 * 1000.0 * 1000.0 * time );

        HOST_PRINT( rank, "Have written " << nBytes << " Bytes" )
        HOST_PRINT( rank, "Write took " << time << " seconds, is " << rateGBs << " GB/s" )
    }
    else
    {
        SCAI_REGION( "main.in" )

        IndexType K;
        IndexType N;

        auto inFile = comm->collectiveFile();

        double time = common::Walltime::get();

        inFile->open( fileName, "r" );
        inFile->readSingle( K );

        HOST_PRINT( rank, "Read collective file " << fileName << " with " << K << " data arrays" )

        hmemo::HArray<ValueType> myData;   // reuse allocated data for reading

        size_t nBytes = 0;

        for ( IndexType iter = 0; iter < K; ++iter )
        {
            SCAI_REGION( "main.read" )

            inFile->readSingle( N );

            auto dist = blockDistribution( N, comm );
            inFile->readAll(  myData, dist->getLocalSize(), dist->lb() );

            // check for correct values in last iteration
            
            if ( iter % 20 == 0 )
            {
                SCAI_REGION( "main.check" )

                auto rData = hostReadAccess( myData );

                for ( IndexType i = 0; i < myData.size(); ++i )
                {
                    SCAI_ASSERT_EQ_ERROR( rData[i], fillArray( dist->local2Global( i  ) ), "illegal result" )
                }
            }

            if ( iter + 1 == K )
            {
                nBytes = inFile->getOffset();
                inFile->close();
                time = common::Walltime::get() - time;
            }
        }

        double rateGBs = nBytes / ( 1000.0 * 1000.0 * 1000.0 * time );

        HOST_PRINT( rank, "Have read " << nBytes << " Bytes" )
        HOST_PRINT( rank, "Read took " << time << " seconds, is " << rateGBs << " GB/s" )
    }
}
