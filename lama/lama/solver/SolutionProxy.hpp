/**
 * @file SolutionProxy.hpp
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
 * @brief Contains the Interface of the Class SolutionProxy
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#ifndef LAMA_SOLUTIONPROXY_HPP_
#define LAMA_SOLUTIONPROXY_HPP_

// for dll_import
#include <common/config.hpp>
#include <common/unique_ptr.hpp>

// others
#include <lama/Vector.hpp>

namespace lama
{

/**
 * @brief The SolutionProxy is used to avoid multiple residual calculations
 *
 * The SolutionProxy is used to avoid needless residual re-calculation.
 * Basically it is a wraps around pointer to a vector. Each not read-only
 * access to the underlying vector marks the proxy as dirty and by that signals
 * the solver to recalculate the residual.
 */
class COMMON_DLL_IMPORTEXPORT SolutionProxy
{
public:

    /**
     * @brief Creates a SolutionProxy with no associated vector.
     */
    SolutionProxy();

    /**
     * @brief Creates a SolutionProxy for the given pointer
     *
     * SolutionProxy does not take owner ship of the passed pointer. Therefore the
     * caller has to make sure that the associated vector exists as long as this
     * SolutionProxy exists.
     *
     * @param[in] solution   The pointer to the vector which the proxy will wrap.
     */
    SolutionProxy( Vector* const solution );

    /**
     * @brief SolutionProxy destructor.
     *
     * The destructor ~SolutionProxy() does not destroy the associated vector.
     */
    ~SolutionProxy();

    /**
     * @brief Returns a constant reference to the underlying vector.
     *
     * @return constant reference to the underlying vector.
     */
    const Vector& getConstReference() const;

    /**
     * @brief Returns a reference to the underlying vector.
     *
     * This call is equivalent to SolutionProxy::getReference(). It was
     * introduced to use the SolutionProxy in the same manner a pointer can
     * be used.
     *
     * @return Reference to the underlying vector.
     */
    Vector& operator*();

    /**
     * @brief Associates the given Vector Pointer with this SolutionProxy.
     *
     * The assignment drops the probably existing association to anther Vector.
     *
     * @param[in] newVector the Vector to which the SolutionProxy shall point to.
     */
    void operator=( Vector* const newVector );

    /**
     * @brief Determines if the proxy is dirty and the residual needs to be
     *        recomputed.
     *
     * @return Value to determine if the proxy is dirty.
     */
    bool isDirty() const;

    /**
     * @brief Sets/unsets the isDirtY flag of the SolutionProxy
     *
     * @param[in] isDirty   value determining whether the proxy is dirty or not.
     */
    void setDirty( bool isDirty );

    /**
     * @brief Returns a reference to the underlying vector.
     *
     * Returns a reference to the underlying vector. One may also use the * operator.
     *
     * @return Reference to the underlying vector.
     */
    Vector& getReference();

    Vector* create();

    void swap( Vector*& other );

private:

    /**
     * @brief The underlying solution vector.
     */
    Vector* mSolution;

    /**
     * @brief Flag which determines, if the Proxy is dirty or not.
     */
    bool mIsDirty;
};

} // namespace lama

#endif // LAMA_SOLUTIONPROXY_HPP_
