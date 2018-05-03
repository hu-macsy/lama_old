/**
 * @file ConvertMatrixStorageBenchmark.hpp
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
 * @brief Benchmark to measure the storage conversion times between different sparse matrix storages.
 * @author Thomas Brandes, Jiri Kraus
 * @date 02.12.2011
 */

#pragma once

#include <scai/benchmark.hpp>
#include <scai/benchmark/Parser.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>


#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace lama
{

/** (Serial) Benchmark to convert matrix storage from one format to another one.
 *
 *  Target and source format of the storage can be specified via the argument
 *  used in the constructor. Value type is given by the template argument.
 *
 *  This benchmark class registers in the benchmark factory via "Convert_ValueType".
 *
 *  @tparam ValueType specifies the type of the matrix storage.
 *
 *  \code
 *     ConvertMatrixStorageBenchmark<float> bench1( "CSR, JDS" );
 *     common::unique_ptr<Benchmark> bench2( Benchmark::create( "Convert_double", "ELL, DIA" ) );
 *     common::unique_ptr<Benchmark> bench3( Benchmark::createWithArgument( "Convert_ComplexFloat", "DENSE, CSR" ) );
 *  \endcode
 */
template<typename ValueType>
class ConvertMatrixStorageBenchmark : 

    public benchmark::Benchmark, 
    private benchmark::Benchmark::Register<ConvertMatrixStorageBenchmark<ValueType> >

{
public:

    ConvertMatrixStorageBenchmark( const std::string& arguments );

    virtual ~ConvertMatrixStorageBenchmark();

    virtual common::ScalarType getValueType() const;

    virtual bool isThreadded() const;

    virtual const std::string& getCreateId() const;

    virtual const std::string& getArgument() const;

    static std::string createValue()
    {
        return std::string( "Convert_" ) + std::string( common::TypeTraits<ValueType>::id() );
    }

    static Benchmark* create( std::string arguments )
    {
        return new ConvertMatrixStorageBenchmark( arguments );
    }

protected:

    virtual void initialize();
    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

private:

    ConvertMatrixStorageBenchmark();

    std::unique_ptr<MatrixStorage<ValueType>> mSourceStorage;
    std::unique_ptr<MatrixStorage<ValueType>> mTargetStorage;

    std::string mArgument;

    SCAI_LOG_DECL_STATIC_LOGGER( logger );
};

}

}
