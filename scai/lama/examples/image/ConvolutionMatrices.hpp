/**
 * @file ConvolutionMatrices.hpp
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
 * @brief Definitions of some convolution matrices to be used for image processing.
 * @author Pascal Maczey
 * @date 13.06.2016
 */
#pragma once

namespace scai
{
namespace imageprocessing
{

const double blur[9] =
{
    0.0, 0.2,  0.0,
    0.2, 0.2,  0.2,
    0.0, 0.2,  0.0
};

const double blur2[25] =
{
    0, 0, 1, 0, 0,
    0, 1, 1, 1, 0,
    1, 1, 1, 1, 1,
    0, 1, 1, 1, 0,
    0, 0, 1, 0, 0,
};

const double findEdges[25] =
{
    0,  0, -1,  0,  0,
    0,  0, -1,  0,  0,
    0,  0,  2,  0,  0,
    0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,
};

const double sharpen[9] =
{
    -1, -1, -1,
    -1,  9, -1,
    -1, -1, -1
};

const double median[9] =
{
    1/9, 1/9, 1/9,
    1/9, 1/9, 1/9,
    1/9, 1/9, 1/9
};

const double sobelX[9] =
{
    1, 0, -1,
    2, 0, -2,
    1, 0, -1
};

const double sobelY[9] =
{
    1,  2,  1,
    0,  0,  0,
    -1, -2, -1
};

const double scharrX[9] =
{
    3, 0,  -3,
    10, 0, -10,
    3, 0,  -3
};

const double scharrY[9] =
{
    3,  10,  3,
    0,   0,  0,
    -3, -10, -3
};

const double gaussian[25] =
{
    0.003765,    0.015019,    0.023792,    0.015019,    0.003765,
    0.015019,    0.059912,    0.094907,    0.059912,    0.015019,
    0.023792,    0.094907,    0.150342,    0.094907,    0.023792,
    0.015019,    0.059912,    0.094907,    0.059912,    0.015019,
    0.003765,    0.015019,    0.023792,    0.015019,    0.003765
};


}  // namespace imageprocessing

}  // namespace scai
