/**
 * @file BenchmarkStub.hpp
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
 * @brief Benchmark Stub that can be used as a copy and paste template to implement new benchmarks.
 * @author Jiri Kraus
 * @date 02.12.2011
 * $Id$
 */
#ifndef LAMA_LAMABENCHMARKSTUB_HPP_
#define LAMA_LAMABENCHMARKSTUB_HPP_

#include <bench/LAMAMPIBenchmark.hpp>

template<typename MatrixType>
class BenchmarkStub: public LAMAMPIBenchmark
{
public:

    typedef typename MatrixType::ValueType ValueType;

    static const std::string& id();

    BenchmarkStub();

    BenchmarkStub( const std::string& arguments );

    BenchmarkStub( const BenchmarkStub<MatrixType>& other );

    virtual ~BenchmarkStub();

    virtual std::auto_ptr<bf::Benchmark> copy() const;

    virtual short getValueTypeSize() const;

    virtual bool isThreadded() const;

    virtual const std::string& getId() const;

protected:
    virtual void initialize();
    virtual void setUp();
    virtual void execute();
    virtual void tearDown();
    virtual void shutdown();

    virtual CounterType getNumFloatingPointOperations() const;
    virtual CounterType getProcessedBytes() const;

    using LAMAMPIBenchmark::mComm;

private:

    static const LAMAInputSetComplexityVisitor::Group& group();

    static const std::string& sid();
};

template<typename MatrixType>
const std::string& BenchmarkStub<MatrixType>::sid()
{
    static const std::string sid = "LAMA<T>";
    return sid;
}

template<typename MatrixType>
const std::string& BenchmarkStub<MatrixType>::id()
{
    static const std::string id = sid() + ": " + LAMAInputSetComplexityVisitor::getGroupId( group() );
    return id;
}

template<typename MatrixType>
BenchmarkStub<MatrixType>::BenchmarkStub()
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ) )
{
}

template<typename MatrixType>
BenchmarkStub<MatrixType>::BenchmarkStub( const std::string& arguments )
    : LAMAMPIBenchmark( sid(), LAMAInputSetComplexityVisitor::getGroupId( group() ), arguments )
{
}

template<typename MatrixType>
BenchmarkStub<MatrixType>::BenchmarkStub( const BenchmarkStub<MatrixType>& other )
    : LAMAMPIBenchmark( other )
{
}

template<typename MatrixType>
BenchmarkStub<MatrixType>::~BenchmarkStub()
{
}

template<typename MatrixType>
std::auto_ptr<bf::Benchmark> BenchmarkStub<MatrixType>::copy() const
{
    bf::Benchmark* b = new BenchmarkStub<MatrixType>( *this );
    return std::auto_ptr<bf::Benchmark>( b );
}

template<typename MatrixType>
short BenchmarkStub<MatrixType>::getValueTypeSize() const
{
    return sizeof(ValueType);
}

template<typename MatrixType>
bool BenchmarkStub<MatrixType>::isThreadded() const
{
    return true;
}

template<typename MatrixType>
const std::string& BenchmarkStub<MatrixType>::getId() const
{
    return id();
}

template<typename MatrixType>
void BenchmarkStub<MatrixType>::initialize()
{
    LAMA_LOG_INFO( logger, "initialize" );

    //device initialization + comunicator

    map < string, string > tokens;

    getConfig (tokens);

}

template<typename MatrixType>
void BenchmarkStub<MatrixType>::setUp()
{
    LAMA_LOG_INFO( logger, "enter Benchmark::setUp" );
}

template<typename MatrixType>
void BenchmarkStub<MatrixType>::execute()
{
    LAMA_LOG_INFO( logger, "execute" );
}

template<typename MatrixType>
void BenchmarkStub<MatrixType>::tearDown()
{
    LAMA_LOG_INFO( logger, "tear down" );
}

template<typename MatrixType>
void BenchmarkStub<MatrixType>::shutdown()
{
    LAMA_LOG_INFO( logger, "shutdown benchmark" );
}

template<typename MatrixType>
CounterType BenchmarkStub<MatrixType>::getNumFloatingPointOperations() const
{
    return 0;
}

template<typename MatrixType>
CounterType BenchmarkStub<MatrixType>::getProcessedBytes() const
{
    return 0;
}

#endif // LAMA_LAMABENCHMARKSTUB_HPP_
