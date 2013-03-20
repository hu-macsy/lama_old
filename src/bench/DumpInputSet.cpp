/**
 * @file DumpInputSet.cpp
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
 * @brief DumpInputSet.cpp
 * @author Jiri Kraus, Robin Rehrmann
 * @date 06.04.2011
 * $Id$
 */

#include <framework/src/benchmark_framework.hpp>

#include <bench/LAMAInputSet.hpp>

extern "C" bf::BaseInputSetRegistry* getInputSetRegistry();

int main( int argc, const char* argv[] )
{
    if( argc != 3 )
    {
        std::cerr << "Wrong number of parameters: " << argv[0] << " <InputSetID> <format>" << std::endl;
        std::cerr << "\t<InputSetID>: The ID of the InputSet to be dumped." << std::endl;
        std::cerr << "\t<format>    : FORMATTED, BINARY, UNFORMATTED, XDR or MATRIX_MARKET" << std::endl;
        return 2;
    }

    lama::File::FileType fileType;
    std::string formatString = argv[2];
    if( formatString == "FORMATTED" )
    {
        fileType = lama::File::FORMATTED;
    }
    else if( formatString == "BINARY" )
    {
        fileType = lama::File::BINARY;
    }
    else if( formatString == "XDR" )
    {
        fileType = lama::File::XDR;
    }
    else if( formatString == "MATRIX_MARKET" )
    {
        fileType = lama::File::MATRIX_MARKET;
    }
    else
    {
        std::cerr << "DumpInputSets: Unknown filetype '" << formatString << '\'' << std::endl;
        return 1;
    }

//  bf::InputSetRegistry<LAMAInputSet>& iSetRegistry = bf::InputSetRegistry<LAMAInputSet>::getRegistry();

    try
    {
        bf::BaseInputSetRegistry& iSetBaseRegistry = ( *getInputSetRegistry() );
        bf::InputSetRegistry<LAMAInputSet>& iSetRegistry =
            dynamic_cast<bf::InputSetRegistry<LAMAInputSet>&>( iSetBaseRegistry );

        const LAMAInputSet& iSet = iSetRegistry.get( argv[1] );
        std::string name = iSet.getId();
        std::cout << "Writing input set " << name << std::endl;

        const lama::CSRSparseMatrix<double>& matrixA = iSet.getA();
        std::cout << " - creating '" << name << '\'' << std::endl;
        matrixA.writeToFile( name, fileType, lama::File::DOUBLE, lama::File::INT, lama::File::INT );

        const lama::DenseVector<double>& vectorX = iSet.getX();
        std::string xname = name + ".sol";
        std::cout << " - creating '" << xname << '\'' << std::endl;
        vectorX.writeToFile( xname, fileType, lama::File::DOUBLE );

        const lama::DenseVector<double>& vectorY = iSet.getY();
        std::string yname = name + ".rhs";
        std::cout << " - creating '" << yname << '\'' << std::endl;
        vectorY.writeToFile( yname, fileType, lama::File::DOUBLE );
    }
    catch( std::exception& e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done!" << std::endl;
    return 0;
}
