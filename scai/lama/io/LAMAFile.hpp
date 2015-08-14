/**
 * @file LAMAFile.hpp
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
 * @brief LAMAFile.hpp
 * @author Thomas Brandes
 * @date 31.10.2011
 * @since 1.0.0
 */
#ifndef LAMA_LAMAFILE_HPP_
#define LAMA_LAMAFILE_HPP_

// for dll_import
#include <scai/common/config.hpp>

namespace lama
{

/** Common base class for LAMA input and output files. */

class COMMON_DLL_IMPORTEXPORT LAMAFile
{
public:

    /**
     * @brief Defines the supported file types
     */
    enum FileType
    {
        /**
         * @brief binary format without header informations in the data file
         */
        BINARY,
        /**
         * @brief ascii format
         */
        FORMATTED,
        /**
         * @brief binary format with header informations in the data file
         */
        UNFORMATTED,
        /**
         * @brief xdr binary format which considers the endianess
         */
        XDR,
        /**
         * @brief the Matrix Market Format
         *        (see http://math.nist.gov/matrixMarket for details).
         */
        MATRIX_MARKET
    };

    LAMAFile( FileType filetype )
        : mFileType( filetype )
    {
    }

protected:

    FileType mFileType;
};

}

#endif // LAMA_LAMAFILE_HPP_
