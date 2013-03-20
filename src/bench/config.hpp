/**
 * @file config.hpp
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
#ifndef LAMABENCH_DLL_IMPORTEXPORT
#ifdef LAMABENCH_COMPILING_DLL
#define LAMABENCH_DLL_IMPORTEXPORT   __declspec( dllexport )
#else //LAMABENCH_COMPILING_DLL is defined
#define LAMABENCH_DLL_IMPORTEXPORT   __declspec( dllimport )
#endif //LAMABENCH_COMPILING_DLL
#endif //LAMABENCH_DLL_IMPORTEXPORT
#else //WIN32 is not defined
#ifndef LAMABENCH_DLL_IMPORTEXPORT
#define LAMABENCH_DLL_IMPORTEXPORT
#endif //LAMABENCH_DLL_IMPORTEXPORT
#endif //WIN32
