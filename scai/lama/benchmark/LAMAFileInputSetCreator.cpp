/**
 * @file LAMAFileInputSetCreator.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief LAMAFileInputSetCreator.cpp
 * @author jiri
 * @date 06.04.2011
 * $Id$
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
