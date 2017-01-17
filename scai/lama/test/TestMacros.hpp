/**
 * @file scai/lama/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

#include <scai/common/test/TestMacros.hpp>

#ifdef SCAI_CHECK_CLOSE
#undef  SCAI_CHECK_CLOSE
#endif

/*
 * @brief HelperMacro SCAI_CHECK_CLOSE( x, y, tolerance )
 *
 * Extends macro BOOST_CHECK_CLOSE( x, y, tolerance) from Boost.Test.
 * Extends macro SCAI_CHECK_CLOSE( x, y, tolerance) from common
 * now works also if x or y is a Scalar
 *
 * @param ValueType     value type to be used for check
 * @param x             Scalar
 * @param y             Scalar
 * @param percent_eps   Epsilon[%]
 *
 */

#define SCAI_CHECK_CLOSE( x, y, tolerance )                                                                   \
    {                                                                                                         \
        scai::lama::Scalar xScalar = scai::lama::Scalar( x );                                                 \
        scai::lama::Scalar yScalar = scai::lama::Scalar( y );                                                 \
        ScalarRepType xVal = xScalar.getValue<ScalarRepType>();                                               \
        ScalarRepType yVal = yScalar.getValue<ScalarRepType>();                                               \
        BOOST_CHECK_CLOSE( scai::common::Math::real( xVal ), scai::common::Math::real( yVal ), tolerance );   \
        BOOST_CHECK_CLOSE( scai::common::Math::imag( xVal ), scai::common::Math::imag( yVal ), tolerance ) ;  \
    }

