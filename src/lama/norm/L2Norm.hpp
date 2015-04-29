/**
 * @file L2Norm.hpp
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
 * @brief L2Norm.hpp
 * @author Jiri Kraus
 * @date 01.06.2011
 * @since 1.0.0
 */
#ifndef LAMA_L2NORM_HPP_
#define LAMA_L2NORM_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/norm/Norm.hpp>

namespace lama
{

/**
 * @brief L2Norm is a functor specialization of Norm that calculates the l2 norm
 *        for the passed values.
 */
class LAMA_DLL_IMPORTEXPORT L2Norm: public lama::Norm
{
public:
    /**
     * @brief Constructs a L2Norm functor.
     */
    L2Norm();

    /**
     * @brief Destroys this L2Norm functor.
     */
    virtual ~L2Norm();

    /**
     * @brief calculates the l2 norm of the passed Scalar.
     *
     * @param[in] scalar    the Scalar to calculate the l2 norm of.
     *
     * @return              the l2 norm of the Scalar scalar.
     */
    virtual Scalar apply( const Scalar& scalar ) const;

    /**
     * @brief calculates the l2 norm of the passed Vector.
     *
     * @param[in] vector    the Vector to calculate the l2 norm of.
     *
     * @return              the l2 norm of the Vector vector.
     */
    virtual Scalar apply( const Vector& vector ) const;

    /**
     * @brief calculates the l2 norm of the passed Matrix.
     *
     * @param[in] matrix    the Matrix to calculate the l2 norm of.
     *
     * @return              the l2 norm of the Matrix matrix.
     */
    virtual Scalar apply( const Matrix& matrix ) const;
};

/**
 * @brief calculates the l2 norm of the passed Scalar.
 *
 * @param[in] scalar    the Scalar to calculate the l2 norm of.
 *
 * @return              the l2 norm of the Scalar scalar.
 */
LAMA_DLL_IMPORTEXPORT Scalar l2Norm( const Scalar& scalar );

/**
 * @brief calculates the l2 norm of the passed Vector.
 *
 * @param[in] vector    the Vector to calculate the l2 norm of.
 *
 * @return              the l2 norm of the Vector vector.
 */
LAMA_DLL_IMPORTEXPORT Scalar l2Norm( const Vector& vector );

/**
 * @brief calculates the l2 norm of the passed Matrix.
 *
 * @param[in] vector    the Matrix to calculate the l2 norm of.
 *
 * @return              the l2 norm of the Matrix matrix.
 */
LAMA_DLL_IMPORTEXPORT Scalar l2Norm( const Matrix& matrix );

}

#endif // LAMA_L2NORM_HPP_
