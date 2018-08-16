/**
 * @file KernelRegistryException.hpp
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
 * @brief Derived exception class to allow specific exceptions thrown at KernelRegistry
 * @author Thomas Brandes
 * @date 10.10.2015
 */

#pragma once

// base class
#include <scai/common/macros/throw.hpp>

namespace scai
{

namespace kregistry
{

/** Derived class needed to catch exception only from KernelRegistry */

class COMMON_DLL_IMPORTEXPORT KernelRegistryException : public scai::common::Exception
{
public:

    KernelRegistryException();

    KernelRegistryException( const std::string& message );

    virtual ~KernelRegistryException() throw();

    virtual const char* what() const throw();

protected:

    std::string mMessage;
};

} /* end namespace common */

} /* end namespace scai */
