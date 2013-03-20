/**
 * @file lama/config.hpp
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
 * @brief config.hpp
 * @author Jiri Kraus
 * @date 20.10.2011
 * $Id$
 */
//No include header guards be cause we want to allow this header to be included multiple times
#ifdef WIN32
#ifdef min
#undef min
#endif //min
#ifdef max
#undef max
#endif //max
//Do not display warnings about dll-interface issues.
// @todo How can we resolve these issues? Do we want to resolve these issues?
#pragma warning( disable: 4251 )

#ifndef LAMA_DLL_IMPORTEXPORT
#ifdef LAMA_COMPILING_DLL
#define LAMA_DLL_IMPORTEXPORT   __declspec( dllexport )
#define LAMA_DLL_EXTERN_TEMPLATE
#else //LAMA_COMPILING_DLL is defined
#define LAMA_DLL_IMPORTEXPORT   __declspec( dllimport )
#define LAMA_DLL_EXTERN_TEMPLATE extern
#endif //LAMA_COMPILING_DLL
#endif //LAMA_DLL_IMPORTEXPORT
#else //WIN32 is not defined
#ifndef LAMA_DLL_IMPORTEXPORT
#define LAMA_DLL_IMPORTEXPORT
#endif //LAMA_DLL_IMPORTEXPORT
#endif //WIN32
