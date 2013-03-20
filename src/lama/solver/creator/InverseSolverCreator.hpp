/**
 * @file InverseSolverCreator.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief InverseSolverCreator.hpp
 * @author Kai Buschulte
 * @date 20.06.2012
 * $Id$
 */
#ifndef LAMA_InverseSolverCreator_HPP_
#define LAMA_InverseSolverCreator_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/creator/SolverCreator.hpp>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT InverseSolverCreator: public SolverCreator
{
public:
    using SolverCreator::RuleType;

    InverseSolverCreator( const std::string type );

    virtual ~InverseSolverCreator();

    /**
     * @brief Returns the unique ID of this solver creator
     */
    static const std::string& id();

    /**
     * @brief Returns the main creator rule of this solver
     */
    virtual RuleType& getCreatorRule();

protected:
    /**
     * @brief Main Rule which is required by the MetaSolver to create a InverseSolver instance
     * MetaSolver is not able to understand the solver-dependent parameters.
     * This is the task of this rule and its sub-rules
     */
    RuleType mRInverseSolver;
private:
    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

} // namespace lama

#endif // LAMA_InverseSolverCreator_HPP_
