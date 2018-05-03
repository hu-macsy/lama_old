/**
 * @file FileInputSet.cpp
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
 * @brief FileInputSet.cpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#include <scai/lama/benchmark/FileInputSet.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/common/Settings.hpp>

namespace scai
{

namespace lama
{

FileInputSet::FileInputSet( const std::string filename ) : LAMAInputSet()
{
    std::string prefix;

    if ( filename.find( "/" ) == 0 )
    { 
        mFileName = filename;
    }
    else if ( common::Settings::getEnvironment( prefix, "SCAI_INPUT_PATH" ) )
    {
        mFileName = prefix + "/" + filename;
    }
    else
    {
        mFileName = filename;
    }

    mAlpha = 1.0;
    mBeta  = 1.0;

    if ( FileIO::fileExists( mFileName ) )
    {
        mA.reset( new CSRSparseMatrix<double>() );
        mA->readFromFile( mFileName );
        mX.reset( new DenseVector<double>() );
        mX->setSameValue( mA->getNumColumns(), 1.0 );
        mY.reset( new DenseVector<double>() );
        *mY = *mA * *mX;    
    }
    else
    {
        SCAI_LOG_ERROR( logger, "file " << mFileName << " does not exist" )

        // do not throw an error, as input set might be created to query arguments

        mFileName = "<valid filename>";
    }
}

benchmark::InputSet* FileInputSet::create( const std::string argument )
{
    // argment is take as correspoding file name

    return new FileInputSet( argument );
}

std::string FileInputSet::createValue()
{
    std::string id = "File";
    return id;
}

const std::string& FileInputSet::getCreateId() const
{
    static std::string id = createValue();
    return id;
}

const std::string& FileInputSet::getArgument() const
{
    return mFileName;
}

}

}
