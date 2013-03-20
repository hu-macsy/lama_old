/**
 * @file Creator.hpp
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
 * @brief Creator.hpp
 * @author Kai Buschulte
 * @date 06.06.2012
 * $Id$
 */
#ifndef LAMA_CREATOR_HPP_
#define LAMA_CREATOR_HPP_

// for dll_import
#include <lama/config.hpp>

#include <logging/Logger.hpp>

#include <lama/Scalar.hpp>

#include <string>

// spirit
#include <boost/spirit/include/qi.hpp>

#include <string>

namespace lama
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

class MetaSolver;

class LAMA_DLL_IMPORTEXPORT Creator
{
public:
    virtual ~Creator();

protected:
    Creator();

    /**
     * @brief Rule to parse a variable name i.e. the ID of a Solver
     */
    qi::rule<std::string::const_iterator,std::string(),ascii::space_type> mRId;

    qi::rule<std::string::const_iterator,Scalar(),ascii::space_type> mRScalar;

private:
    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

}
#endif // LAMA_CREATOR_HPP_
