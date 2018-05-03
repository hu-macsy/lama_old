/**
 * @file utilskernel/examples/vectorAdd.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./utilskernel/examples/vectorAdd.cpp
 * @author Lauretta Schubert
 * @date 29.02.2016
 */
#include <scai/hmemo.hpp>
#include <scai/utilskernel.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

using namespace scai;

void add( hmemo::HArray<double> res, const hmemo::HArray<double> a, const hmemo::HArray<double> b )
{
    SCAI_ASSERT_LE( a.size(), b.size(), "size mismatch" )
    IndexType n = a.size();
    hmemo::ContextPtr hostCtx = hmemo::Context::getContextPtr( common::ContextType::Host );
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

int main( int, char** )
{
    int size = 10;
    hmemo::HArray<double> a, b, c;
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<double> > setVal;
    hmemo::ContextPtr loc = hmemo::Context::getContextPtr( common::ContextType::Host );
    setVal.getSupportedContext( loc );
    hmemo::WriteOnlyAccess<double> writeB( b, loc, size );
    hmemo::WriteOnlyAccess<double> writeC( c, loc, size );
    setVal[loc]( writeB.get(), size, double( 2 ), common::BinaryOp::COPY );
    setVal[loc]( writeC.get(), size, double( 3 ), common::BinaryOp::COPY );
    writeB.release();
    writeC.release();
    add ( a, b, c );
    exit ( 0 );
}
