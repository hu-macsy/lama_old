/**
 * @file CUDASettings.hpp
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
 * @brief Managing some settings for CUDA specified by environment variables
 * @author Thomas Brandes
 * @date 04.05.2013
 * @since 1.0.0
 */

#pragma once

// internal scai libraries
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** This class determines whether the texture and or shared memory on a GPU
 *  device should be used or not.
 *
 *  This class provides only static methods.
 */

class CUDASettings
{

public:

    /** Check whether textures should be used or not in kernel algorithms. */

    static bool useTexture();

    /** Check wheter shared memory should be used or not in kernel algorithms. */

    static bool useSharedMem();

    /** Get the block size, might be set by environment variable. */

    static int getBlockSize();

    /** Get the block size, n is the degree of parallelism */

    static int getBlockSize( const int n );

    /** Enable or disable texture and shared memory use explicitly. */

    static void set( bool useSharedMemFlag, bool useTextureFlag );

private:

    CUDASettings();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized; //!< will be set true after determination of theUseTextureFlag

    static bool theUseTextureFlag;//!< result for useTexture if initialized is true

    static bool theUseSharedMemFlag;//!< result for useTexture if initialized is true

    static int theBlockSize;//!< result for blockSize if initialized is true

    /**
     *   Get the (major) compute capability for the current active device.
     */
    static int getComputeCapability();

    /**
     *   Initialization of flags.
     */
    static void initialize();
};

} /* end namespace lama */

} /* end namespace scai */
