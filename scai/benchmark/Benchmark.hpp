/**
 * @file Benchmark.hpp
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
 * @brief Definition of a common base class for all benchmarks.
 * @author Jiri Kraus, Thomas Brandes
 * @date 04.05.2010, revised 14.09.2017
 */

#pragma once

#include <string>
#include <list>
#include <algorithm>
#include <numeric>
#include <memory>

#include <scai/common/config.hpp>
#include <scai/common/Factory1.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Printable.hpp>

#include <scai/benchmark/BenchmarkTypes.hpp>

#include <scai/logging.hpp>

namespace scai
{

/**
 * @brief The namespace benchmark holds everything of the benchmark framework
 */
namespace benchmark
{

/** Common base class for all benchmark classes.
 *
 *  New benchmark classes must derive from this base class and should register themselves
 *  at the benchmark factory.
 *
 *  Each derived benchmark class must implement certain routines so that it can be
 *  called in a unified manner.
 */
class COMMON_DLL_IMPORTEXPORT Benchmark : 

    public common::Factory1< std::string, std::string, Benchmark* >,
    public common::Printable,
    private common::NonCopyable

{
public:

    /**
     * @brief Constructor of the base class
     *
     * @param[in] name      name of the benchmark used in statistics
     * @param[in] groupName group name of the benchmark, used for grouping results
     */
    Benchmark( const std::string& name, const std::string& groupId );

    /**
     * @brief The destructor destroys this object and frees all inner resources.
     */
    virtual ~Benchmark();

    /** Override default implementation of Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;
    /**
     * @brief executes the Benchmark
     *
     * This function runs the Benchmark by calling its protected functions
     * initialize( ), setUp( ), execute( ) (getNumRepititions( )-times),
     * tearDown( ) and shutdown( ), whereas the time of executing of each
     * setUp( ), execute( ) and tearDown( ) is saved.
     *
     * @param[in] out   The outputstream, used for writing out information of
     *                  running the Benchmark.
     */
    virtual void run( std::ostream& out );

    /**
     * @brief Returns the name of the benchmark.
     * @return The name of the benchmark.
     */
    virtual const std::string& getName() const;

    /**
     * @brief Returns the name of the group, the benchmark belongs to.
     * @return The name of the group.
     */
    virtual const std::string& getGroup() const;

    /**
     * @brief Returns the create id of the benchmark, i.e. the class name of the derived benchmark class.
     *
     * @return id that has been used to create the benchmark.
     */
    virtual const std::string& getCreateId() const = 0;

    /**
     * @brief Return the argument string that has been used to create the benchmark.
     */
    virtual const std::string& getArgument() const = 0;

    /**
     * @brief Returns the number of threads, if threadded, 'unthreadded' otherwise.
     * @return The number of threads.
     */
    virtual std::string getNumThreads() const;

    /**
     * @brief Returns information, whether the platform is adapted for multi
     *        threadding.
     * @return true, if this platform supports multi threadding, false,
     * otherwise.
     */
    virtual bool isThreadded() const;

    /**
     * @brief Sets the time, the benchmark will be executed minimal.
     *
     * This function sets the minimum time for the execution of the function
     * execute( ). If the function execute( ) is run and the given time is not
     * over, the execution will be repeated - this is true for the number of
     * repititions, also. However, a running execution will not be interrupted,
     * that is, why it is the minimal time, not the absolute.
     *
     * @param[in] minTime - the time the function execute( ) will be executed,
     *                      at least. Unit is seconds.
     */
    virtual void setMinTime( const double minTime );

    /**
     * @brief Returns the number of repititions.
     * @return The number of repititions.
     */
    int getNumRepitions() const;

    /**
     * @brief Sets the number of repititions.
     *
     * The number of repititions defines, how many times the function execute( )
     * is called.
     *
     * @param[in] numRepitions The number of repititions.
     */
    void setNumRepitions( const int numRepitions );

    /**
     * @brief Sets the id of the InputSet of this benchmark.
     * @param[in] inputSetId The id of the InputSet.
     */
    virtual void setInputSetId( const std::string& inputSetId );

    /**
     * @brief Returns the id of the InputSet.
     * @return The id of the InputSet.
     */
    virtual const std::string& getInputSetId() const;

    /**
     * @brief Returns the setup time.
     * @return The time it took, to execute setup( ).
     */
    virtual double getSetupTime() const;

    /**
     * @brief Returns the mean average execution time.
     *
     * Computes the mean average of all times it took to run execute( ) in
     * consideration of outlier.
     *
     * @return The mean average time it took, to run execute( ).
     */
    virtual double getExecutionTime() const;

    /**
     * @brief Returns the minimal time it took, to run execute( ).
     * @return The minimal time it took, to run execute( ). At least 0.0.
     */
    virtual double getMinExecutionTime() const;

    /**
     * @brief Returns the maximal time it took, to run execute( ).
     * @return The minimal time it took, to run execute( ). 0.0 if execute( )
     * was
     *         not run.
     */
    virtual double getMaxExecutionTime() const;

    /**
     * @brief Returns the number of floating points per second. The unit of the
     *        result is giga.
     * @return The number of 1.000.000.000 floating points per second.
     */
    virtual double getExecutionFlops() const;

    /**
     * @brief Returns the number of bytes processed per second. The unit of the
     *        result is giga.
     * @return The number of 1.073.741.824 bytes processed per second.
     */
    virtual double getExecutionBandwidth() const;

    /**
     * @brief Returns the teardown time.
     * @return The time it took, to execute teardown( ).
     */
    virtual double getTearDownTime() const;

    /**
     * @brief Returns the benchmark execution time.
     * @return The time it took, to execute the benchmark.
     */
    virtual double getTotalExecutionTime() const;

    /**
     * @brief Returns the value type of benchmark data.
     */
    virtual common::ScalarType getValueType() const = 0;

    /**
     * @brief Returns the time, function execute( ) was run.
     * @return the time the function execute( ) was run in seconds.
     */
    virtual double getExecutedTime() const;

    /**
     * @brief Returns the number of repititions of executing execute( ).
     * @return The number of repititions of executing execute( ).
     */
    virtual int getNumActualRepititons() const;

    /**
     * @brief This struct is used to sort benchmarks according to their group
     *        ids.
     */
    struct CompareGroupIds
    {
        bool operator( )( const Benchmark* const arg1, const Benchmark* const arg2 );
    };
    /**
     * @brief This struct is used to sort benchmarks according to their
     *        expression ids.
     */
    struct CompareInputSetIds
    {
        bool operator( )( const Benchmark* const arg1, const Benchmark* const arg2 );
    };

    /**
     * @brief This struct is used to sort benchmarks according to their
     *        expression ids.
     */
    struct CompareValueTypeSize
    {
        bool operator( )( const Benchmark* const arg1, const Benchmark* const arg2 );
    };

    /**
     * @brief This struct is used to sort benchmarks according to their
     *        number of threads.
     */
    struct GreaterNumThreads
    {
        bool operator( )( const Benchmark* const arg1, const Benchmark* const arg2 );
    };

    /**
     * @brief This struct is used to sort benchmarks according to their
     *        number of threads.
     */
    struct LessNumThreads
    {
        bool operator( )( const Benchmark* const arg1, const Benchmark* const arg2 );
    };

    virtual bool doOutput() const;

    /** @brief create of a input set via "Benchmark( argument )" */

    static Benchmark* createWithArgument( const std::string& specification );

protected:

    /**
     * @brief Gets additional information by the given string.
     *
     * A Benchmark may has additional flags to be set or different options,
     * which would be usually given by comandlineparameters. These comandline-
     * parameters are passed with this string argument. It leaves the user the
     * opportunity to add information to the benchmark, by parsing the given
     * string.
     * The default behaviour, though, is, to throw an exception, because the
     * abstract class does not support this feature.
     * This function is called by the two constructors, offering the option to
     * pass arguments with a string. If the user wishes, to use this feature,
     * all that has to be done, is, to overwrite this function.
     *
     * @param[in] arguments The string, holding the arguments for the benchmark.
     * @throw               Exception if this feature is not supported by this
     *                      benchmark.
     */
    virtual void interpreteArgs( const std::string& arguments );

    /**
     * @brief Creates preconditions for the preparation and execution of the
     *        benchmark.
     *
     * Creates preconditions for the preparation and execution of the
     * benchmark. For example this function could read benchmark-data from an
     * external file.
     * The execution of this function is not timed. Therefore it can prepare
     * anything for the benchmark that does not have to be timed.
     */
    virtual void initialize() = 0;

    /**
     * @brief Prepares the benchmark.
     *
     * Prepares the benchmark. For example this function could initialize
     * important variables or upload data from host to device.
     */
    virtual void setUp() = 0;

    /**
     * @brief Executes the main part of the benchmark.
     *
     * Executes the main part of the benchmark. For example this function could
     * execute a complicate algorithm.
     */
    virtual void execute() = 0;

    /**
     * @brief Cleans memory.
     *
     * Cleans memory. For example this function could free those variables,
     * initialized in setUp( ) or download data from device to host.
     */
    virtual void tearDown() = 0;

    /**
     * @brief Cleans up everything, before leaving the program.
     *
     * Cleans up everything, before leaving the program. For example this
     * function could close open streams to files or anywhere else.
     * The execution of this function is not timed. Therefore it can clean
     * anything from the benchmark that does not have to be timed.
     */
    virtual void shutdown() = 0;

    virtual void synchronize() const;

    /**
     * @brief Returns the absolute number of floating point operations.
     *
     * @todo: getNumFloatingPointOperations can have a default implementation
     *        so we should add it.
     * @return The absolute number of floating point operations.
     */
    virtual CounterType getNumFloatingPointOperations() const = 0;

    /**
     * @brief Returns the absolute number of processed bytes.
     *
     * @todo: getProcessedBytes can have a default implementation
     *        so we should add it.
     * @return The absolute number of processed bytes.
     */
    virtual CounterType getProcessedBytes() const = 0;

    /** Name of this benchmark, used in statistics */
    std::string mName;  

    std::string mGroup;  //!< The name of the group

    /** The number of threads on which the benchmark will be executed. */
    int mNumThreads;
    /** The minimum runtime. */
    double mMinTime;
    /** The time of execute( ) being run. */
    double mExecutedTime;
    /** The id of the InputSet. */
    std::string mInputSetId;
    /** The final number of repititions. */
    int mActNumRep;

    /** Logger for Benchmark, can be used in derived classes. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

private:

    /**
     * @brief Disable the default constructor.
     */
    Benchmark();

    int mNumRepitions;
    double mSetUpTime;
    double mTearDownTime;
    std::list<double> mExecutionTimes;

};

} //namespace benchmark

}
