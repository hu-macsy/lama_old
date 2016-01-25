/**
 * @file SOR.hpp
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
 * @brief SOR.hpp
 * @author Martin Schenk
 * @date 29.08.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lamasolver/Solver.hpp>
#include <scai/lamasolver/OmegaSolver.hpp>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT SOR:
		public OmegaSolver,
		public Solver::Register<SOR>
{
public:
    SOR( const std::string& id );

    SOR( const std::string& id, const Scalar omega );

    SOR( const std::string& id, LoggerPtr logger );

    SOR( const std::string& id, const Scalar omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    SOR( const SOR& other );

    virtual ~SOR();

    void initialize( const Matrix& coefficients );

    void iterate();

    struct SORRuntime: OmegaSolverRuntime
    {
        SORRuntime();
        virtual ~SORRuntime();
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual SORRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SORRuntime& getConstRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    SORRuntime mSORRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

private:
    template<typename ValueType>
    void iterateImpl();

    common::unique_ptr<const Matrix> mIterationMatrix;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
