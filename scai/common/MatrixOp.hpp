/**
 * @file MatrixOp.hpp
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
 * @brief Enum type for some common matrix operations.
 * @author Thomas Brandes
 * @date 27.02.2018
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/macros/throw.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
{

/** Enumeration type for matrix operators used implicitly in matrix-vector multiplication
 */
enum class MatrixOp
{
    NORMAL,          //!< no operations at all
    CONJ,            //!< uses conj of all matrix elements
    TRANSPOSE,       //!< use tranpose of the matrix
    CONJ_TRANSPOSE,  //!< use conjugate-transpose of the matrix
    MAX_MATRIX_OP    //!< internal use only
};

inline bool isTranspose( const MatrixOp op )
{
    return op == MatrixOp::TRANSPOSE || op == MatrixOp::CONJ_TRANSPOSE;
}

inline bool isConj( const MatrixOp op )
{
    return op == MatrixOp::CONJ || op == MatrixOp::CONJ_TRANSPOSE;
}

/*
 * Output of MatrixOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const MatrixOp& op )
{
    switch ( op )
    {
        case MatrixOp::NORMAL:
            stream << "NORMAL";
            break;

        case MatrixOp::CONJ:
            stream << "CONJ";
            break;

        case MatrixOp::TRANSPOSE:
            stream << "TRANSPOSE";
            break;

        case MatrixOp::CONJ_TRANSPOSE:
            stream << "CONJ_TRANSPOSE";
            break;

        default:
            stream << "<unknown_matrix_op " << static_cast<int>( op ) << ">";
            break;
    }

    return stream;
}

/** Combination of two matrix operations 
 *
 *  @param[in] op1 first operation
 *  @param[in] op2 second operation
 *  @returns   the combination of the two operations
 *
 *  \code
 *    matrixA = transpose( conj( B ) )
 *  \endcode
 */
inline MatrixOp combine( const MatrixOp op1, const MatrixOp op2 )
{
    MatrixOp op = MatrixOp::MAX_MATRIX_OP;

    if ( op1 == MatrixOp::NORMAL )
    {
        op = op2;
    }
    else if ( op2 == MatrixOp::NORMAL )
    {
        op = op1;
    }
    else if ( op1 == op2 )
    {
        // same operation applied twice is identity 

        op = MatrixOp::NORMAL;
    }
    else if ( MatrixOp::CONJ == op1 )
    {
        if ( MatrixOp::TRANSPOSE == op2 )
        {  
            op = MatrixOp::CONJ_TRANSPOSE;
        }
        else if ( MatrixOp::CONJ_TRANSPOSE == op2 )
        {
            op = MatrixOp::TRANSPOSE;
        }
    }
    else if ( MatrixOp::TRANSPOSE == op1 )
    {
        if ( MatrixOp::CONJ == op2 )
        {  
            op = MatrixOp::CONJ_TRANSPOSE;
        }
        else if ( MatrixOp::CONJ_TRANSPOSE == op2 )
        {
            op = MatrixOp::CONJ;
        }
    }
    else if ( MatrixOp::CONJ_TRANSPOSE == op1 )
    {
        if ( MatrixOp::CONJ == op2 )
        {  
            op = MatrixOp::TRANSPOSE;
        }
        else if ( MatrixOp::TRANSPOSE == op2 )
        {
            op = MatrixOp::CONJ;
        }
    }

    if ( op == MatrixOp::MAX_MATRIX_OP )
    {
        COMMON_THROWEXCEPTION( "combine( " << op1 << ", " << op2 << " undefined" )
    }

    return op;
}

} /* end namespace common */

} /* end namespace scai */
