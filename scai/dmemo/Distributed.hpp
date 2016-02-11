/**
 * @file Distributed.hpp
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
 * @brief Distributed.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/dmemo/Distribution.hpp>

namespace scai
{

namespace dmemo
{

/** Common base class for all objects that are distributed.
 *
 * Important: There is no default constructor, a distribution must
 * always be specified. You can use NoDistribtion for non distributed objects.
 * */

class COMMON_DLL_IMPORTEXPORT Distributed: public scai::common::Printable
{
public:

    Distributed( DistributionPtr );

    Distributed( const Distributed& other );

    virtual ~Distributed();

    inline const Distribution& getDistribution() const;

    inline DistributionPtr getDistributionPtr() const;

    /** This method returns connectivity information that will be used to repartition with Metis.
     *
     *  @param[out] ia offsets for ja array, size is localSize+1
     *  @param[out] ja connectivities, for elem i it is ja[ia[i]], ..., ja[ia[i+1]-1]
     *  @param[out] vwgt are weights, size is localSize
     *
     *  Note: the size of the ja array is given get getCSRGraphSize()
     */

    virtual void buildCSRGraph( IndexType ia[], IndexType ja[], IndexType vwgt[], const IndexType* globalRowIndexes ) const;

    /** This method returns the number of connectivities, becomes size of ja array in buildCSRGraph */

    virtual IndexType getCSRGraphSize() const;

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

} /* end namespace dmemo */

} /* end namespace scai */
