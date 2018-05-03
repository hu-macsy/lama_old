/**
 * @file HArrayBenchmark.hpp
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
 * @brief Benchmark for copying HArray data between different devices
 * @author Thomas Brandes, Jiri Kraus
 * @date 19.09.2017
 */

#pragma once

#include <scai/benchmark.hpp>

#include <scai/hmemo/benchmark/HArrayInputSet.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>

namespace scai
{

namespace hmemo
{

/** Benchmark to measure transfer time between two devices.
 *
 *  This benchmark registers at the benchmark factory, benchmarks can be created
 *  by HArrayBenchmark( targetContext, sourceContext ).
 *
 *  This benchmark must be used with an HArrayInputSet as input set.
 *
 *  \code
 *     HArrayBenchmark inputData( "Host, CUDA" );
 *     std::unique_ptr<HArrayBenchmark> inputData1( Benchmark::create( "HArrayBenchmark", "CUDA, Host" ) );
 *     std::unique_ptr<HArrayBenchmark> inputData2( Benchmark::createWithArg( "HArrayBenchmark( Host, CUDA )" ) );
 *  \endcode
 */
class HArrayBenchmark: 
 
    public  benchmark::Benchmark,
    private benchmark::Benchmark::Register<HArrayBenchmark>

{
public:

    /** Constructor of the benchmark.
     *
     *  @param[in] argument must be a string like "CSR, double"
     */
    HArrayBenchmark( const std::string& argument );

    virtual ~HArrayBenchmark();

    /** Implementation of pure method Benchmark::getCreateId()   */

    virtual const std::string& getCreateId() const
    {
        static std::string id = createValue();
        return id;
    }

    /** Implementation of pure method Benchmark::getArgument()   */

    virtual const std::string& getArgument() const
    {
        return mArgument;
    }

    static std::string createValue()
    {
        return "HArrayBenchmark";
    }

    static Benchmark* create( const std::string argument )
    {
        return new HArrayBenchmark( argument );
    }

    virtual common::ScalarType getValueType() const
    {
        return common::ScalarType::DOUBLE;
    }

    /** Override implementation Benchmark::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

    HArrayBenchmark();

    /** Implementation of pure method Benchmark::initialize */

    virtual void initialize();

    /** Implementation of pure method Benchmark::setUp */

    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

    std::string mArgument;  // argument used for creating this benchmark

    std::unique_ptr<benchmark::InputSet> mInputSet;

    HArrayInputSet* mHArrayInputSet;

    HArray<double>* mArray;   // pointer to the input array
    
    ContextPtr mTargetContext;
    ContextPtr mSourceContext;

};

}

}
