/**
 * @file SBLASRandomInputSetCreator.cpp
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
 * @brief LAMARandomInputSetCreator.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */
/*
 * LAMARandomInputSetCreator.cpp
 *
 *  Created on: 17.11.2010
 *      Author: dfeld
 */

#include <bench/LAMARandomInputSetCreator.hpp>
#include <sstream>
#include <cstdlib>

#include <iostream>

using lama::IndexType;

namespace
{

static double randomNumber()
{
    return rand() / static_cast<double>( RAND_MAX );
}

}

const std::string& LAMARandomInputSetCreator::id()
{
    static const std::string id =
        bf::InputSetRegistry<LAMARandomInputSetCreator::InputSetType>::getRandomInputSetCreatorId();
    return id;
}

LAMARandomInputSetCreator::LAMARandomInputSetCreator()
{
}

LAMARandomInputSetCreator::~LAMARandomInputSetCreator()
{
}

std::auto_ptr<LAMARandomInputSetCreator::InputSetType> LAMARandomInputSetCreator::create() const
{
    throw bf::BFException( "LAMARandomInputSetCreator creator needs parameters to create a InputSet." );
}

std::auto_ptr<LAMARandomInputSetCreator::InputSetType> LAMARandomInputSetCreator::create(
    const std::string& params ) const
{
    srand( static_cast<unsigned int>( time( NULL ) ) );

    std::istringstream input( params );

    IndexType size = 0;
    double fillingGrade = 1.0;

    input >> size;

    if( size <= 0 ) {
        throw bf::BFException( "LAMARandomInputSetCreator creator needs a size to create a InputSet." );
    }

    input >> fillingGrade;

    if( input.good() && ( fillingGrade <= 0.0 || fillingGrade > 1.0 ) )
        throw bf::BFException(
            "LAMARandomInputSetCreator creator needs valid filling grade > 0.0 and <= 1.0 to create a InputSet." );

    if( input.fail() )
    {
        fillingGrade = 1.0;
    }

    std::auto_ptr<lama::DenseMatrix<double> > matrixA( new lama::DenseMatrix<double>( size, size ) );

    double value = 1.0 / static_cast<double>( size );

    lama::HostWriteAccess<double> matrixValues( matrixA->getLocalStorage().getData() );

    if( fillingGrade < 1.0 )
    {
        #pragma omp parallel for schedule(static)
        for( IndexType i = 0; i < size; ++i )
        {
            for( IndexType j = 0; j < size; ++j )
            {
                if( randomNumber() <= fillingGrade )
                {
                    matrixValues[i * size + j] = value;
                }
                else
                {
                    matrixValues[i * size + j] = 0.0;
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for schedule(static)
        for( IndexType i = 0; i < size; ++i )
        {
            for( IndexType j = 0; j < size; ++j )
            {
                matrixValues[i * size + j] = value;
            }
        }
    }

    std::auto_ptr<lama::DenseVector<double> > vectorX( new lama::DenseVector<double>( size, 1.0 ) );
    std::auto_ptr<lama::DenseVector<double> > vectorY( new lama::DenseVector<double>( size, 1.0 ) );

    lama::HostWriteAccess<double> y( vectorY->getLocalValues() );
    lama::HostReadAccess<double> x( vectorX->getLocalValues() );

    if( fillingGrade < 1.0 )
    {
        #pragma omp parallel for schedule(static)
        for( IndexType i = 0; i < size; ++i )
        {
            y[i] = 0.0;
            for( IndexType j = 0; j < size; ++j )
            {
                y[i] += matrixValues[i * size + j] * x[j];
            }
        }
    }

    std::auto_ptr<InputSetType> iSet( new InputSetType( getId(), 1.0, 1.0, vectorX, vectorY, matrixA ) );

    return iSet;
}

const std::string& LAMARandomInputSetCreator::getId() const
{
    return id();
}

LAMA_INPUTSET_REGISTRATION( LAMARandomInputSetCreator );
