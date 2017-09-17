/**
 * @file FileInputSet.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief FileInputSet.hpp
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#pragma once

#include <scai/lama/benchmark/LAMAInputSet.hpp>
#include <scai/benchmark.hpp>

#include <string>

namespace scai
{

namespace lama
{

/** Derived class of LAMAInputSet that creates the data by an input file */

class FileInputSet: public  LAMAInputSet, 
                    private benchmark::InputSet::Register<FileInputSet>
{
public:

    /** Constructor of a file input set
     *
     *  This method does not throw an error if the filename does not exist but the
     *  LAMAInputSet remains empty.
     */
    FileInputSet( const std::string filename );

    /** Implemenation of method that returns the key for creating an object via the factory. */

    static std::string createValue();

    /** Implementation of pure method InputSet::getCreateId()   */

    virtual const std::string& getCreateId() const;

    /** Static method required for create, called by InputSet::create( InputSet::createValue(), argument ) */

    static InputSet* create( const std::string argument );

    /**
     * @brief Implementation of pure method InputSet::getArgument()
     */
    virtual const std::string& getArgument() const;

private:

    std::string mFileName;
};

}

}
