/**
 * @file Distributed.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
