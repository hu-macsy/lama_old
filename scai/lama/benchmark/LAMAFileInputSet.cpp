/**
 * @file LAMAFileInputSet.cpp
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
 * @brief LAMAFileInputSet.cpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#include <scai/lama/benchmark/LAMAFileInputSet.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

namespace scai
{

namespace lama
{

LAMAFileInputSet::LAMAFileInputSet( const std::string filename ) : LAMAInputSet()
{
    std::string file;

    if ( filename.find( "/" ) == 0 )
    { 
        file = filename;
    }
    else
    {
        std::string prefix = benchmark::Config::getInstance().getValueOf( "path" );
        file = prefix + filename;
    }

    mAlpha = 1.0;
    mBeta  = 1.0;

    mA.reset( new CSRSparseMatrix<double>( file ) );
    mX.reset( new DenseVector<double>( mA->getNumColumns(), 1.0 ) );
    mY.reset( new DenseVector<double>( *mA * *mX  ) );
}

benchmark::InputSet* LAMAFileInputSet::create( const std::string argument )
{
    // argment is take as correspoding file name

    return new LAMAFileInputSet( argument );
}

std::string LAMAFileInputSet::createValue()
{
    std::string id = "File";
    return id;
}

}

}
