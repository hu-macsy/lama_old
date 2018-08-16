/**
 * @file CUDASettings.hpp
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
 * @brief Managing some settings for CUDA specified by environment variables
 * @author Thomas Brandes
 * @date 04.05.2013
 */

#pragma once

#include <scai/common/config.hpp>

namespace scai
{

namespace common
{

/** This class determines whether the texture and or shared memory on a GPU
 *  device should be used or not.
 *
 *  This class provides only static methods.
 */

class COMMON_DLL_IMPORTEXPORT CUDASettings
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

} /* end namespace common */

} /* end namespace scai */
