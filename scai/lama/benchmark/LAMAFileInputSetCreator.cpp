/**
 * @file LAMAFileInputSetCreator.cpp
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
 * @brief LAMAFileInputSetCreator.cpp
 * @author jiri
 * @date 06.04.2011
 */
/**
 * @file LAMAFileInputSetCreator.cpp
 * @author jiri
 * Created on: 10.08.2010
 */

#include <scai/lama/benchmark/LAMAFileInputSetCreator.hpp>
#include <scai/lama/benchmark/LAMAInputSetComplexityVisitor.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

using namespace scai;

const std::string& LAMAFileInputSetCreator::id()
{
    static const std::string id =
        bf::InputSetRegistry<LAMAFileInputSetCreator::InputSetType>::getFileInputSetCreatorId();
    return id;
}

LAMAFileInputSetCreator::LAMAFileInputSetCreator()
{
}

LAMAFileInputSetCreator::~LAMAFileInputSetCreator()
{
}

LAMAFileInputSetCreator::InputSetType* LAMAFileInputSetCreator::create() const
{
    throw bf::BFError( "File LAMAFileInputSetCreator creator needs a filename to create a InputSet." );
}

LAMAFileInputSetCreator::InputSetType* LAMAFileInputSetCreator::create( const std::string& filename ) const
{
    std::string prefix = bf::Config::getInstance().getValueOf( "path" );

    std::string file = prefix + filename;

    common::unique_ptr<lama::CSRSparseMatrix<double> > sparseMatrix( new lama::CSRSparseMatrix<double>( file ) );

    common::unique_ptr<lama::DenseVector<double> > VectorX(
        new lama::DenseVector<double>( sparseMatrix->getNumColumns(), 1.0 ) );

    common::unique_ptr<lama::DenseVector<double> > VectorY(
        new lama::DenseVector<double>( ( *sparseMatrix ) * ( *VectorX ) ) );

    return new LAMAInputSet( getId(), 1.0, 1.0, VectorX, VectorY, sparseMatrix );
}

const std::string& LAMAFileInputSetCreator::getId() const
{
    return id();
}

LAMA_INPUTSET_REGISTRATION( LAMAFileInputSetCreator );
