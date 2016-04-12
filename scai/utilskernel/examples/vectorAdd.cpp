#include <scai/hmemo.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;

void add( hmemo::HArray<double> res, const hmemo::HArray<double> a, const hmemo::HArray<double> b )
{
    SCAI_ASSERT_LE( a.size(), b.size(), "size mismatch" )

    IndexType n = a.size();
    hmemo::ContextPtr hostCtx = hmemo::Context::getContextPtr( common::context::Host );

    hmemo::ReadAccess<double> read1( a, hostCtx ); 
    hmemo::ReadAccess<double> read2( b, hostCtx );
    hmemo::WriteOnlyAccess<double> write( res, hostCtx, n );

    double* resPtr = write.get();
    const double* aPtr = read1.get();
    const double* bPtr = read2.get();

    for ( IndexType i = 0; i < n; ++i )
    {
        resPtr[i] = aPtr[i] + bPtr[i];
    }
}

int main(int, char**)
{
    int size = 10;
    hmemo::HArray<double> a, b, c;

    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<double, double> > setVal;

    hmemo::ContextPtr loc = hmemo::Context::getContextPtr( common::context::Host );

    setVal.getSupportedContext( loc );

    hmemo::WriteOnlyAccess<double> writeB( b, loc, size);
    hmemo::WriteOnlyAccess<double> writeC( c, loc, size);

    setVal[loc]( writeB.get(), size, double( 2 ), common::reduction::COPY );
    setVal[loc]( writeC.get(), size, double( 3 ), common::reduction::COPY );

    writeB.release();
    writeC.release();

    add ( a, b, c );

    exit (0);
}
