/**
 * @file Distributed.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Distributed.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
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

    // TODO: document globalRowIndexes
    /** This method returns connectivity information that will be used to repartition with Metis.
     *
     *  @param[out] ia               offsets for ja array, size is localSize+1
     *  @param[out] ja               connectivities, for elem i it is ja[ia[i]], ..., ja[ia[i+1]-1]
     *  @param[out] vwgt             are weights, size is localSize
     *  @param[out] globalRowIndexes 
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
