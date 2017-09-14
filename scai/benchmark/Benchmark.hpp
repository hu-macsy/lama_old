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
 * @brief Benchmark.hpp
 * @author jiri
 * @date 04.05.2010
 */

#pragma once

#include <string>
#include <list>
#include <algorithm>
#include <numeric>
#include <memory>
#include <omp.h>

#include <scai/common/config.hpp>
#include <scai/benchmark/BenchmarkTypes.hpp>

#include <scai/logging.hpp>

/**
 * @brief The namespace bf holds everything of the benchmark framework
 */
namespace bf
{
class COMMON_DLL_IMPORTEXPORT Benchmark
{
public:
    /**
     * @brief The constructor creates a default Benchmark-object with default
     *        values.
     */
    Benchmark();

    /** Overrride the default copy constructor as timing is reset */
    Benchmark( const Benchmark& other );

    /**
     * @brief The constructor creates a Benchmark-object.
     * @param[in] id    The id of this Benchmark.
     * @param[in] gid   The id of the group of this Benchmark.
     */
    Benchmark( const std::string& id, const std::string& gid );

private:

    // ToDo: seems to be unused, arguments is ignored

    /**
     * @brief creates a Benchmark with the given arguments.
     * @param[in] arguments A string of arguments for the benchmark.
     */
    Benchmark( const std::string& arguments );

public:

    /**
     * @brief The constructor creates a Benchmark-object.
     * @param[in] id        The id of this Benchmark.
     * @param[in] gid       The id of the group of this Benchmark.
     * @param[in] arguments A string of arguments for the benchmark.
     */
    Benchmark( const std::string& id, const std::string& gid, const std::string& arguments );

    /**
     * @brief The destructor destroys this object and frees all inner resources.
     */
    virtual ~Benchmark();

    /**
     * @brief Copies this Benchmark.
     * @return The copy of this benchmark.
     */
    virtual Benchmark* copy() const = 0;

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
    virtual const std::string& getGid() const;

    /**
     * @brief Returns the ID of the benchmark.
     * @return The id of the benchmark.
     */
    virtual const std::string& getId() const = 0;

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
     * @param[in] numRepititions The number of repititions.
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
     * @brief Returns the size of the value type.
     * @return The size of the value type.
     */
    virtual short getValueTypeSize() const = 0;

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
     * @brief Compares this benchmark to another.
     *
     * Two benchmarks are the same, if they have the same ID, because the ID of
     * the benchmark makes it unique.
     *
     * @param[in] other The other benchmark.
     * @return true, if both benchmarks are identical, false, otherwise.
     */
    virtual bool operator==( const Benchmark& other );

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

    class HasId
    {
    public:
        HasId( const std::string& id );
        bool operator( )( const Benchmark* const benchmark );
    private:
        std::string mId;
    };

    virtual bool doOutput() const;

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

    /** The id of the expression. */
    std::string mName;
    /** The id of the group. */
    std::string mGId;
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
    int mNumRepitions;
    double mSetUpTime;
    double mTearDownTime;
    std::list<double> mExecutionTimes;

};

} //namespace bf
