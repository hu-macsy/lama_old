/**
 * @file PMVBenchmark.hpp
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
 * @brief Benchmark for parallel matrix-vector multiplication with different settings.
 * @author Thomas Brandes, Jiri Kraus
 * @date 19.09.2017
 */

#pragma once

#include <scai/lama/benchmark/LAMAMPIBenchmark.hpp>
#include <scai/lama/benchmark/LAMAInputSet.hpp>

namespace scai
{

namespace lama
{

/** (Parallel) benchmark for matrix-vector multiplication.
 *
 *  This benchmark registers at the benchmark factory, benchmarks can be created
 *  by PMV( <matrix format>, <value type> ).
 *
 *  This benchmark can be used with any LAMAInputSet, i.e. with all objects that
 *  derive from this class.
 */
class PMVBenchmark: 
 
    public  LAMAMPIBenchmark,
    private benchmark::Benchmark::Register<PMVBenchmark>

{
public:

    typedef double ValueType;     // should become a template argument

    /** Constructor of the benchmark.
     *
     *  @param[in] argument must be a string for the format like "CSR"
     */
    PMVBenchmark( const std::string& argument );

    virtual ~PMVBenchmark();

    virtual common::ScalarType getValueType() const;

    virtual bool isThreadded() const;

    /** Implementation of pure method Benchmark::getCreateId()   */

    virtual const std::string& getCreateId() const;

    /** Implementation of pure method Benchmark::getArgument()   */

    virtual const std::string& getArgument() const;

    static std::string createValue();

    static Benchmark* create( const std::string argument )
    {
        return new PMVBenchmark( argument );
    }

protected:

    /** Implementation of pure method Benchmark::initialize */

    virtual void initialize();

    /** Implementation of pure method Benchmark::setUp */

    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

    using benchmark::Benchmark::mName;

    using LAMAMPIBenchmark::mComm;

private:

    PMVBenchmark();

    static const std::string& getGroupId();

    static void getComplexity( 
        CounterType& numFlops, 
        CounterType& numProcessedIndexes, 
        CounterType& numPprocessedValues,
        const Matrix<ValueType>& matrix );

    std::unique_ptr<benchmark::InputSet> mInputSet;
    const LAMAInputSet* mLAMAInputSet;

    std::unique_ptr<Matrix<ValueType>> mMatrixA;
    std::unique_ptr<DenseVector<ValueType>> mVectorX;
    std::unique_ptr<DenseVector<ValueType>> mVectorY;

    hmemo::ContextPtr mContext;
    SyncKind mCommunicationKind;

    common::ScalarType mType;  // value type of input data

    CounterType mNumFloatingPointOperations;
    CounterType mNumProcessedIndexes;
    CounterType mNumProcessedValues;

    std::string mArgument;  // argument used for creating this benchmark
};

}

}
