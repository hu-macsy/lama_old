/**
 * @file LAMAMPIBenchmark.hpp
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
 * @brief LAMAMPIBenchmark.hpp
 * @author Jiri Kraus
 * @date 13.05.2011
 * $Id$
 */
#ifndef LAMA_LAMAMPIBENCHMARK_HPP_
#define LAMA_LAMAMPIBENCHMARK_HPP_

#include <framework/src/Benchmark.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <map>
#include <string>

class LAMAMPIBenchmark: public bf::Benchmark
{
public:

    LAMAMPIBenchmark();

    LAMAMPIBenchmark( const std::string& id, const std::string& gid );

    LAMAMPIBenchmark( const std::string& id, const std::string& gid, const std::string& arguments );

    LAMAMPIBenchmark( const LAMAMPIBenchmark& other );

    virtual ~LAMAMPIBenchmark();

    LAMAMPIBenchmark& operator=( const LAMAMPIBenchmark& other );

    virtual std::string getNumThreads() const;

    virtual bool doOutput() const;

    /** This method delivers for each MPI process a vector of
     *  configuration tokens.
     *
     *  The configuration is set by the environment variable LAMA_CONFIG.
     *  The tokens are separated by a colon.
     *
     *  \code
     *    mpirun -np 4 -x LAMA_CONFIG="HOST=4:W=3"
     *  \endcode
     *
     *  It is also possible to specify a different configuration for each
     *  process.
     *
     *  \code
     *    mpirun -np 3 -x LAMA_CONFIG="HOST=4:W=3,CUDA=0:W=6,CUDA=1:W=6"
     *  \endcode
     *
     *  Configurations for different processors are separated by a comma.
     *  If more than one configuration is specified, the number must match
     *  the number of processors.
     *
     *  \@param[out] tokens is map of config string key value pairs,
     *               e.g. <"LOCAL"="HOST", "THREADS"="4", "W"="3">
     *
     *  The interpretation of tokens is completely left to the calling routine.
     */
    void getConfig( std::map<std::string,std::string>& tokens ) const;

protected:

    virtual void synchronize() const;

    lama::CommunicatorPtr mComm;

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

private:

    void readConfig();

    std::string config;
};

#endif // LAMA_LAMAMPIBENCHMARK_HPP_
