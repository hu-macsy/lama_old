/**
 * @file lama_info.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Example program that prints internal information about LAMA factories
 * @author Thomas Brandes
 * @date 25.07.2015
 */

// Define levels for assertion, logging and tracing

#include "lama.hpp"

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/expression/all.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matutils/MatrixCreator.hpp>
#include <lama/NoCommunicator.hpp>

#include <common/shared_ptr.hpp>

#include <iostream>

using namespace std;

void contextInfo()
{
    using namespace memory;

    vector<ContextType> values;  // supported context types

    Context::getCreateValues( values );

    cout << "Factory of Context: " << values.size() << " entries" << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "Registered values[" << i << "] = " << values[i] << endl;

        // ContextFactory has additional argument for device id

        ContextPtr context = Context::create( values[i], -1 );
      
        cout << "Context is " << *context << endl;

        cout << "Memory is " << *context->getMemoryPtr() << endl;
        cout << "Host Memory is " << *context->getHostMemoryPtr() << endl;
    }
}

void communicatorInfo()
{
    using namespace lama;

    CommunicatorPtr nocomm = NoCommunicator::create();

    cout << "NoCommunicator : " << *nocomm << endl;

    vector<string> values;  // string is create type for the factory

    Communicator::getCreateValues( values );

    cout << "Factory of Communicator: " << values.size() << " entries" << endl;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "   Registered values[" << i << "] = " << values[i] << endl;
    }

    CommunicatorPtr comm = Communicator::get();
}

void matrixInfo()
{
    using namespace lama;

    vector<MatrixCreateKeyType> keys;

    Matrix::getCreateValues( keys );

    cout << "Factory of Matrix: " << keys.size() << " entries" << endl;

    for ( size_t i = 0; i < keys.size(); ++i )
    {
        cout << "   Registered values[" << i << "] = " << keys[i].first << ", " << keys[i].second << endl;

        common::shared_ptr<Matrix> matrix ( Matrix::create( keys[i] ) );

        cout << "Matrix : " << *matrix << endl;
    }
}

void vectorInfo()
{
    using namespace lama;

    vector<VectorCreateKeyType> keys;

    Vector::getCreateValues( keys );

    cout << "Factory of Vector: " << keys.size() << " entries" << endl;

    for ( size_t i = 0; i < keys.size(); ++i )
    {
        cout << "   Registered values[" << i << "] = " << keys[i].first << ", " << keys[i].second << endl;

        common::shared_ptr<Vector> vector ( Vector::create( keys[i] ) );

        cout << "Vector : " << *vector << endl;
    }
}

int main( int argc, char* argv[] )
{
    communicatorInfo();
    contextInfo();
    vectorInfo();
    matrixInfo();
}

