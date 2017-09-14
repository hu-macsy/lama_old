/**
 * @file readDir.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief readDir.h
 * @author Jiri Kraus
 * @date 06.04.2011
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

namespace scai
{

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

} // namespace scai
