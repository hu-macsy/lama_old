/**
 * @file VTInterface.hpp
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
 * @brief Interface between LAMA and VampirTrace
 * @author Thomas Brandes
 * @date 14.02.2013
 * $Id$
 */
#ifndef LAMA_VT_INTERFACE_HPP_
#define LAMA_VT_INTERFACE_HPP_

namespace tracing
{

struct RegionEntry;

/** Static class provding interface to Vampir */

class VTInterface
{
public:

    /** Type definition needed as region entry gets additional entry to save region id given by VampirTrace. */

    typedef unsigned int VTRegionId;

    /** This routine defines a new region within VampirTrace. */

    static void define( RegionEntry& );

    static void enter( const RegionEntry& region );

    static void leave( const RegionEntry& region );

    static void enable( bool flag );
};

} // namespace

#endif // LAMA_VT_INTERFACE_HPP
