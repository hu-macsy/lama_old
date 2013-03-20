/**
 * @file SORCreator.hpp
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
 * @brief SORCreator.hpp
 * @author Kai Buschulte
 * @date 20.06.2012
 * $Id$
 */
#ifndef LAMA_SORCreator_HPP_
#define LAMA_SORCreator_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/creator/OmegaSolverCreator.hpp>

namespace lama
{

namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class LAMA_DLL_IMPORTEXPORT SORCreator: public OmegaSolverCreator
{
public:
    using OmegaSolverCreator::RuleType;

    SORCreator( const std::string type );

    virtual ~SORCreator();

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
     * @brief Main Rule which is required by the MetaSolver to create a SOR instance
     * MetaSolver is not able to understand the solver-dependent parameters.
     * This is the task of this rule and its sub-rules
     */
    RuleType mRSOR;

private:
    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

} // namespace lama

#endif // LAMA_SORCreator_HPP_
