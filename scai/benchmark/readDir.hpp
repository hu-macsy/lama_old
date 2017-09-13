/**
 * @file readDir.h
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
 * @brief readDir.h
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */
/*
 * readDir.h
 *
 *  Created on: 31.01.2011
 *      Author: rrehrman
 */

#pragma once

#include <vector>
#include <string>

namespace bf
{

/**
 * @brief Returns the names of the shared libraries listed in path in the given
 *        vector.
 * @param[in]  path     The path holding the shared libraries.
 * @param[out] files    The vector holding the names of the shared libraries
 *                      within the given path.
 * @throws Exception    If directory within the path could not be opened.
 */
void getFilesFromPath( const std::string& path, std::vector<std::string>& files );

} // namespace bf

