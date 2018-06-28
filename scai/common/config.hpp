/**
 * @file common/config.hpp
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
 * @brief General configuration file to deal with different features of Windows/Linux
 * @author Jiri Kraus
 * @date 10.06.2015
 */

//No include header guards because we want to allow this header to be included multiple times

#ifdef WIN32

#ifdef min
#undef min
#endif //min

#ifdef max
#undef max
#endif //max

//Do not display warnings about dll-interface issues.
//TODO: How can we resolve these issues? Do we want to resolve these issues?
#pragma warning( disable : 4251 )

#ifndef COMMON_DLL_IMPORTEXPORT
#ifdef COMMON_COMPILING_DLL
#define COMMON_DLL_IMPORTEXPORT   __declspec( dllexport )
#else //COMMON_COMPILING_DLL is defined
#define COMMON_DLL_IMPORTEXPORT   __declspec( dllimport )
#endif //COMMON_COMPILING_DLL
#endif //COMMON_DLL_IMPORTEXPORT

#else  //WIN32

// visibility can now be restricted in gnu compilers with -fvisibility=hidden

#ifndef COMMON_DLL_IMPORTEXPORT

#if __GNUC__ >= 4
#define COMMON_DLL_IMPORTEXPORT  __attribute__ ( ( visibility ( "default" ) ) )
#else
#define COMMON_DLL_IMPORTEXPORT
#endif

#endif


#endif //WIN32

