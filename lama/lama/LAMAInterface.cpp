/**
 * @file LAMAInterface.cpp
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
 * @brief Default implementations for LAMA interface.
 * @author Thomas Brandes
 * @date 28.04.2011
 * @since 1.0.0
 */

// hpp
#include <lama/LAMAInterface.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( LAMAInterface::logger, "LAMAInterface" )

LAMAInterface::LAMAInterface()
{
    // default constructors do all NULL initialization of
    // function pointer variables
}

LAMAInterface::~LAMAInterface()
{
}

UtilsInterface::UtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

CSRUtilsInterface::CSRUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

ELLUtilsInterface::ELLUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

JDSUtilsInterface::JDSUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

DIAUtilsInterface::DIAUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

COOUtilsInterface::COOUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

DenseUtilsInterface::DenseUtilsInterface()
{
    memset( this, 0, sizeof( *this ) );
}

void LAMAInterface::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "LAMAInterface(base class)";
}

} //namespace lama
