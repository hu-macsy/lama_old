/**
 * @file HaloBuilder.hpp
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
 * @brief HaloBuilder.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 24.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_HALOBUILDER_HPP_
#define LAMA_HALOBUILDER_HPP_

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/distribution/Halo.hpp>
#include <scai/lama/distribution/Distribution.hpp>

namespace scai
{

namespace lama
{

class COMMON_DLL_IMPORTEXPORT HaloBuilder
{
public:
    static void build( const Distribution& distribution, const std::vector<IndexType>& requiredIndexes, Halo& halo );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */

#endif // LAMA_HALOBUILDER_HPP_
