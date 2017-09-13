/**
 * @file LAMAFileInputSetCreator.hpp
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
 * @brief LAMAFileInputSetCreator.hpp
 * @author jiri
 * @date 06.04.2011
 * $Id$
 */

#pragma once

#include <scai/benchmark.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>

#include <string>

class LAMAFileInputSetCreator: public bf::InputSetCreator<LAMAInputSet>
{
public:

    typedef bf::InputSetCreator<LAMAInputSet>::InputSetType InputSetType;

    static const std::string& id();

    LAMAFileInputSetCreator();

    virtual ~LAMAFileInputSetCreator();

    virtual InputSetType* create() const;

    virtual InputSetType* create( const std::string& filename ) const;

    virtual const std::string& getId() const;
};
