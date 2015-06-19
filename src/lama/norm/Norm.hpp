/**
 * @file Norm.hpp
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
 * @brief Norm.hpp
 * @author Jiri Kraus
 * @date 01.06.2011
 * @since 1.0.0
 */
#ifndef LAMA_NORM_HPP_
#define LAMA_NORM_HPP_

// for dll_import
#include <common/config.hpp>

// others
#include <lama/matrix/Matrix.hpp>
#include <lama/Vector.hpp>
#include <lama/Scalar.hpp>

namespace lama
{

class Norm;

typedef boost::shared_ptr<Norm> NormPtr;

/**
 * @brief Norm is a abstract functor to calculate the norms for the passed values.
 *
 * The Functor Norm is mainly used by the stopping criteria ResidualThreshold and
 * ResidualStagnation, to allow a to customize the Norms for these stopping criteria.
 */
class COMMON_DLL_IMPORTEXPORT Norm
{
public:
    /**
     * @brief Constructs a empty norm object.
     */
    Norm();

    /**
     * @brief Destroys this.
     */
    virtual ~Norm();

    /**
     * @brief Calculates the norm of the passed Scalar.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] scalar    the Scalar to caluclate the norm for.
     * @return              the norm of scalar.
     */
    Scalar operator()( const Scalar& scalar ) const;

    /**
     * @brief Calculates the norm of the passed Vector.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] vector    the Vector to caluclate the norm for.
     * @return              the norm of vector.
     */
    Scalar operator()( const Vector& vector ) const;

	/**
     * @brief Calculates the norm of the passed Matrix.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] vector    the Matrix to caluclate the norm for.
     * @return              the norm of matrix.
     */
    Scalar operator()( const Matrix& matrix ) const;

    /**
     * @brief Calculates the norm of the passed Scalar.
     *
     * @param[in] scalar    the Scalar to caluclate the norm for.
     * @return              the norm of scalar.
     */
    virtual Scalar apply( const Scalar& scalar ) const =0;

    /**
     * @brief Calculates the norm of the passed Vector.
     *
     * @param[in] vector    the Vector to caluclate the norm for.
     * @return              the norm of vector.
     */
    virtual Scalar apply( const Vector& vector ) const =0;

    /**
     * @brief Calculates the norm of the passed Matrix.
     *
     * @param[in] vector    the Matrix to caluclate the norm for.
     * @return              the norm of matrix.
     */
    virtual Scalar apply( const Matrix& matrix ) const =0;
};

}

#endif // LAMA_NORM_HPP_
