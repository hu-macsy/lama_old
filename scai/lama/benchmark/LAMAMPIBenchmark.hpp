/**
 * @file LAMAMPIBenchmark.hpp
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
 * @brief LAMAMPIBenchmark.hpp
 * @author Jiri Kraus
 * @date 13.05.2011
 */
#pragma once

#include <scai/benchmark/Benchmark.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <map>
#include <string>

namespace scai
{

class LAMAMPIBenchmark: public benchmark::Benchmark
{
public:

    LAMAMPIBenchmark( const std::string& name, const std::string& gId );

    virtual ~LAMAMPIBenchmark();

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
    void getConfig( std::map<std::string, std::string>& tokens ) const;

protected:

    virtual void synchronize() const;

    scai::dmemo::CommunicatorPtr mComm;

    SCAI_LOG_DECL_STATIC_LOGGER( logger );

private:

    void readConfig();

    std::string mConfig;
};

}
