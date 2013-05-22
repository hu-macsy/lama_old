/**
 * @file Distributed.hpp
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
 * @brief Distributed.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_DISTRIBUTED_HPP_
#define LAMA_DISTRIBUTED_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/distribution/Distribution.hpp>

namespace lama
{

/** Common base class for all objects that are distributed.
 *
 * Important: There is no default constructor, a distribution must
 * always be specified. You can use NoDistribtion for non distributed objects.
 * */

class LAMA_DLL_IMPORTEXPORT Distributed: public Printable
{
public:

    Distributed( DistributionPtr );

    Distributed( const Distributed& other );

    virtual ~Distributed();

    inline const Distribution& getDistribution() const;

    inline DistributionPtr getDistributionPtr() const;

protected:

    void swap( Distributed& other );

    void setDistributionPtr( DistributionPtr distributionPtr );

private:

    DistributionPtr mDistribution; // distribution of obj, never NULL

    Distributed(); // disable default constructor

};

const Distribution& Distributed::getDistribution() const
{
    return *mDistribution;
}

DistributionPtr Distributed::getDistributionPtr() const
{
    return mDistribution;
}

}

#endif // LAMA_DISTRIBUTED_HPP_
