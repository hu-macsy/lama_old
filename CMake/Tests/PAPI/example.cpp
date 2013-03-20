#define NUMEVENTS 2

#include <papi.h>
#include <iostream>
#include <sstream>
#include <exception>

#define PAPI_CALL_EXPECTED( call, retval_expected )  \
    {   \
        int retval = call;  \
        if ( retval != retval_expected )  \
        {  \
            std::ostringstream errorStr;                                                \
            errorStr << "PAPI error in line " << __LINE__;                              \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call;                                              \
            errorStr << "  Error: ";                                                    \
            errorStr << PAPI_strerror( retval );                                        \
            errorStr << ", retval = " << retval << std::endl;                           \
            std::cerr <<  errorStr.str();                                               \
            throw std::exception();                                                     \
       }  \
   }
    

#define PAPI_CALL( call ) PAPI_CALL_EXPECTED( call, PAPI_OK )

void compute( int n )
{
    float* a = new float[n];
    float* b = new float[n];
    float* c = new float[n];

    for ( int i = 0; i < n; i++ ) 
    {
        a[i] = 1.0;
        b[i] = 2.5;
    }
    for ( int i = 0; i < n; i++ ) 
    {
        c[i] = a[i] + b[i];
    }
    float s = 0;
    for ( int i = 0; i < n; i++ ) 
    {
        s = s + c[i];
    }

    delete [] a;
    delete [] b;
    delete [] c;

    std::cout << "Result = " << s << std::endl;
}

int main ( int argc, char* argv[] )
{
    int events[NUMEVENTS] = { PAPI_TOT_INS, PAPI_TOT_CYC };
    int errorcode;
    long long values[NUMEVENTS];

    PAPI_CALL_EXPECTED( PAPI_library_init( PAPI_VER_CURRENT ), PAPI_VER_CURRENT );

    PAPI_CALL( PAPI_start_counters( events, NUMEVENTS ) );

    PAPI_CALL( PAPI_read_counters( values, NUMEVENTS ) );

    std::cout << "Counters = " << values[0] << ", " << values[1] << std::endl;

    compute( 1000 );

    PAPI_CALL( PAPI_read_counters( values, NUMEVENTS ) );

    std::cout << "Counters = " << values[0] << ", " << values[1] << std::endl;

    compute( 1000 );

    PAPI_CALL( PAPI_read_counters( values, NUMEVENTS ) );

    std::cout << "Counters = " << values[0] << ", " << values[1] << std::endl;

    PAPI_CALL ( PAPI_stop_counters( values, NUMEVENTS ) );

    std::cout << "Counters = " << values[0] << ", " << values[1] << std::endl;
}
