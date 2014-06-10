/**
 * @file NoDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief NoDistribution.hpp
 * @author brandes
 * @date 14.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_NODISTRIBUTION_HPP_
#define LAMA_NODISTRIBUTION_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/distribution/Distribution.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** Distribution class that stands for a replicated distribution.
 *
 *  With this distribution an object has an incarnation on each 
 *  processor. 
 *
 *  Usually, methods should take care of consistency among 
 *  all processors, i.e. writes and update operations must be
 *  done on all partitions. But a replicated object can also be used
 *  like a private incarnation on each processor.
 */
class LAMA_DLL_IMPORTEXPORT NoDistribution: public Distribution
{
public:

    NoDistribution( const IndexType globalSize );

    virtual ~NoDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual IndexType getLocalSize() const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    void printDistributionVector( std::string name ) const;

    /** Static methods to create a NoDistribution. */

    static NoDistribution* create( const CommunicatorPtr commPtr, const IndexType globalSize, const float weight = 1.0 );

    static NoDistribution* create( const CommunicatorPtr commPtr, const Matrix& matrix, const float weight = 1.0 );

private:

    NoDistribution(); // no default constructor as global size is not available

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static bool initialized;  //!< used for static registration
};

}

#endif // LAMA_NODISTRIBUTION_HPP_
