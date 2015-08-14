/**
 * @file GeneralDistribution.hpp
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
 * @brief GeneralDistribution.hpp
 * @author brandes
 * @date 25.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_GENERALDISTRIBUTION_HPP_
#define LAMA_GENERALDISTRIBUTION_HPP_

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/LAMATypes.hpp>

#include <scai/lama/distribution/Distribution.hpp>

#include <vector>

#if (BOOST_VERSION < 103600)
#include <map>
#else //(BOOST_VERSION >= 103600)
#include <boost/unordered_map.hpp>
#endif //(BOOST_VERSION < 103600)
namespace lama
{

/** A general distribution allows to map a global range of values
 to the partitions of a communicator completely arbitrarily.

 Each partition has the information which values it holds.
 A partition has no information where to find values
 not owned by itself.
 */

class COMMON_DLL_IMPORTEXPORT GeneralDistribution: public Distribution
{
public:

    /** Construcor of a general distribution.
     *  \param globalSize is the size of the distributed range
     \param myGlobalIndexes contains all indexes of range owned by this partition
     \param communicator partitions on which the range is distributed.

     Important: each global index from 0 to globalSize-1 must appear exactly once in
     the vector myGlobalIndexes on one partition.
     */

    GeneralDistribution(
        const IndexType globalSize,
        const std::vector<IndexType>& myGlobalIndexes,
        const CommunicatorPtr communicator );

    GeneralDistribution(
        const std::vector<IndexType>& row2Partition,
        const IndexType globalSize,
        const CommunicatorPtr communicator );

    explicit GeneralDistribution( const Distribution& other );

//    GeneralDistribution(const GeneralDistribution& other);

    virtual ~GeneralDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual IndexType getLocalSize() const;

    virtual std::vector<IndexType>& getLocalRows();

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    void getDistributionVector( std::vector<IndexType>& row2Partition ) const;

    void printDistributionVector( std::string name ) const;

protected:

    GeneralDistribution( const IndexType globalSize, const CommunicatorPtr communicator );

    //TODO: Evaluate which one is faster
#if (BOOST_VERSION < 103600)
    typedef std::map<IndexType,IndexType> Global2LocalMapType;
#else //(BOOST_VERSION >= 103600)
    typedef boost::unordered_map<IndexType,IndexType> Global2LocalMapType;
#endif //(BOOST_VERSION < 103600)
    Global2LocalMapType mGlobal2Local;
    std::vector<IndexType> mLocal2Global;

private:

    GeneralDistribution();

    GeneralDistribution& operator=( const GeneralDistribution& other );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

typedef common::shared_ptr<GeneralDistribution> GeneralDistributionPtr;

}

#endif // LAMA_GENERALDISTRIBUTION_HPP_
