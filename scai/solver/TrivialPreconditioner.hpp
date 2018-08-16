/**
 * @file TrivialPreconditioner.hpp
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
 * @brief TrivialPreconditioner.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT TrivialPreconditioner:

    public Solver<ValueType>

{
public:

    TrivialPreconditioner( const std::string& id );
    TrivialPreconditioner( const std::string& id, LoggerPtr logger );
    TrivialPreconditioner( const TrivialPreconditioner& other );

    virtual ~TrivialPreconditioner();

    virtual void solveImpl();

    virtual TrivialPreconditioner* copy();

    struct TrivialPreconditionerRuntime: Solver<ValueType>::SolverRuntime
    {
        TrivialPreconditionerRuntime();
        virtual ~TrivialPreconditionerRuntime();
    };

    virtual TrivialPreconditionerRuntime& getRuntime();

    virtual const TrivialPreconditionerRuntime& getRuntime() const;

    static std::string createValue();

    static Solver<ValueType>* create( const std::string name );

protected:
    TrivialPreconditionerRuntime mTrivialPreconditionerRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;
};

} /* end namespace solver */

} /* end namespace scai */
