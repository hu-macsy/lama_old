/**
 * @file Lib/config.hpp
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
 * @brief Test of dllimport/dllexport for Windows Libs.
 * @author Thomas Brandes
 * @date 01.05.2013
 */

#ifdef DLL_LINKING

    // Microsoft Visual C++ compiler requires dllimport / dllexport

    #ifdef COMPILING_DLL

        #define DLL_IMPORTEXPORT   __declspec( dllexport )

    #else // COMPILING_DLL is defined

        #define DLL_IMPORTEXPORT   __declspec( dllimport )

    #endif //COMPILING_DLL

#else

    // ignore DLL_IMPORTEXPORT for other compilers

    #define DLL_IMPORTEXPORT 

#endif

