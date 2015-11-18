
#include <scai/kregistry/KernelContextFunction.hpp>

#include <iostream>

using namespace scai;
using namespace scai::kregistry;
using namespace scai::common::context;

struct UtilsInterface
{
    template<typename ValueType>
    struct isSorted
    {
        typedef bool ( *FuncType ) ( const ValueType array[], const int n, bool ascending );
        static inline const char* getId() { return "Utils.IsSorted"; }
    };
};

static bool isSorted( const double a[], int N, bool ascending )
{
    std::cout << "isSorted<double>, N = " << N << ", ascending = " << ascending << std::endl;

    bool is = true;
    for ( int i = 0; i < N-1; ++i )
    {
        if ( ascending )
        {
            is = a[i] <= a[i+1];
        }
        else
        {
            is = a[i] >= a[i+1];
        }
        if ( !is )
        {
            break;
        }
    }
 
    return is;
}

template<typename ValueType>
void scale( ValueType mValues[], const ValueType value, const int n )
{
    if( value == static_cast<ValueType>( 1 ) )
    {
        return;
    }

    if ( value == 0 )
    {
        for( int i = 0; i < n; i++ )
        {
            mValues[i] = 0;
        }
    }
    else
    {
        for( int i = 0; i < n; i++ )
        {
            mValues[i] *= value;
        }
    }
}

static void setInterface()
{
    std::cout << std::endl;
    std::cout << "setInterface: start" << std::endl;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;

    KernelRegistry::set<UtilsInterface::isSorted<double> >( isSorted, Host, flag );
    KernelRegistry::set( isSorted, "Utils.isSorted", Host, flag );

    KernelRegistry::set( scale<float>, "Utils.scale", Host, flag );
    KernelRegistry::set( scale<float>, "Utils.scale", CUDA, flag );
    KernelRegistry::set( scale<double>, "Utils.scale", Host, flag );
    KernelRegistry::set( scale<double>, "Utils.scale", CUDA, flag );

    std::cout << "setInterface: done" << std::endl;
    KernelRegistry::printAll();
}

static void example1()
{
    std::cout << std::endl;
    std::cout << "Example 1:" << std::endl;
    std::cout << "==========" << std::endl;

    typedef bool ( *SigIsSorted ) ( const double*, int N, bool ascending );

    SigIsSorted isSorted;

    // double ( *isSorted ) ( const double*, int, bool );
    
    KernelRegistry::get( isSorted, "Utils.isSorted", Host );

    double a[] = { 3.0, 4.0, 5.0 };
  
    bool okay = isSorted( a, 3, true );

    std::cout << "example1: isSorted = " << okay << std::endl;
}

static void example2()
{
    std::cout << std::endl;
    std::cout << "Example 2:" << std::endl;
    std::cout << "==========" << std::endl;

    /* Alternative use ( not recommended for LAMA ):

        typedef bool ( *SigIsSorted ) ( const double*, int, bool );
    
        static KernelContextFunction< SigIsSorted > isSorted( "Utils.isSorted" );
    */

    KernelTraitContextFunction< UtilsInterface::isSorted<double> > isSorted;

    double a[] = { 3.0, 4.0, 2.0 };
  
    bool okay = isSorted[ Host ]( a, 3, true );

    std::cout << "example2: isSorted = " << okay << std::endl;
}

template<typename ValueType>
static void example3()
{
    std::cout << std::endl;
    std::cout << "Example 3:" << std::endl;
    std::cout << "==========" << std::endl;

    std::cout << "example3, ValueType = " << typeid( ValueType ).name() << std::endl;

    typedef void ( *SigScale ) ( ValueType*, ValueType, int );

    ValueType a[] = { 3, 4, 2 };
    
    static KernelContextFunction< SigScale > scale ( "Utils.scale" ) ;

    scale[ Host ]( a, 10, 3 );

    std::cout << "example3: scale: " << a[0] << ", " << a[1] << ", " << a[2] << std::endl;
}

static void example4()
{
    std::cout << std::endl;
    std::cout << "Example 4:" << std::endl;
    std::cout << "==========" << std::endl;

    static KernelContextFunction< bool (*) ( const double*, int, bool ) > isSorted( "Utils.isSorted" );
    static KernelContextFunction< void (*) ( double*, double, int ) > scale( "Utils.scale" );

    std::cout << "isSorted: valid context = " << isSorted.validContext( CUDA ) << std::endl;
    std::cout << "scale: valid context = " << scale.validContext( CUDA ) << std::endl;
    std::cout << "scale, isSorted: valid context = " << scale.validContext( isSorted, CUDA ) << std::endl;

    std::cout << std::endl;
}

int main()
{
    setInterface();

    example1();
    example2();

    example3<float>();
    example3<double>();

    example4();
}
