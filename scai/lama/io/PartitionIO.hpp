/**
 * @file PartitionIO.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief IO support for partitioned read/write of vectors, matrices, distributions
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#pragma once

#include <scai/dmemo/Distribution.hpp>
#include <scai/logging.hpp>

#include <cstring>

namespace scai
{

namespace lama
{

/** This class contains some routines to support parallel IO by read and write
 *  of distributed data into partitions, i.e. each processor reads and writes its
 *  local data instead of one single processor reads and writes all global data.
 */

class COMMON_DLL_IMPORTEXPORT PartitionIO
{
public:

    static dmemo::DistributionPtr readDistribution(
        const std::string& inFileName,
        dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr() );

    /** This method reads a distribution from a single input file.
     *
     *  This input file must contain an array, where the size is the global size of the array
     *  and the elements are owner for each index.
     *
     *  @param[in] inFileName is the name of the input file containing the owners for each element
     *  @param[in] comm is the Communicator for which the distribution is determined.
     *  @returns a new distribution for the mapping specified in the file
     *
     *  The global size of the distribution is given by the number of entries in the file.
     *  All owners must be values between 0 and size-1 where size is the size of comm.
     *
     *  If the owners are sorted, a general block distribution is returned, otherwise a
     *  general distribution.
     */
    static dmemo::DistributionPtr readSDistribution( const std::string& inFileName, dmemo::CommunicatorPtr comm );

    /** This method reads a distribution from one input file for each partition.
     *
     *  @param[in] inFileName is the name of the input file containing the indexes for this partition
     *  @param[in] comm is the Communicator for which the distribution is determined.
     *  @returns a new distribution for the mapping specified in the file
     *
     *  The global size of the distribution is given by summing up the number of entries in all files.
     *  Every index between 0 and global size - 1 must appear exactly once in one of the files.
     */
    static dmemo::DistributionPtr readPDistribution( const std::string& inFileName, dmemo::CommunicatorPtr comm );

    static void write( const dmemo::Distribution& distribution, const std::string& fileName );

    static void writePDistribution( const dmemo::Distribution& distribution, const std::string& inFileName );

    static void writeSDistribution( const dmemo::Distribution& distribution, const std::string& inFileName );

    /** Help routine to get the filename of the partition that is read/written by this processor.
     *
     *  @param[in,out] fileName is the name of the file that might contain %r that will be substituted by rank_size
     *  @param[out] isPartitioned if true data is spread among multiple 'partitoned' files
     *  @param[in] comm is the communicator needed to get the rank if name contains %r
     */

    static void getPartitionFileName( std::string& fileName, bool& isPartitioned, const dmemo::Communicator& comm );

    static void getPartitionFileName( std::string& fileName, bool& isPartitioned, const PartitionId rank, const PartitionId size );

    static void getSingleFileName( std::string& fileName );

    /** Counterpart to FileIO::removeFile with partitioned file names.
     *
     *  @param[in] fileName file to delete
     *  @param[in] comm is the communicator
     *  @returns   0 on success, same result on all processors
     */
    static int removeFile( const std::string& fileName, const dmemo::Communicator& comm );

    static bool isPartitionFileName( const std::string& fileName );

    /** Counterpart to FileIO::fileExists, now with partitioned file names */

    static bool fileExists( const std::string& fileName, const dmemo::Communicator& comm );

    static const PartitionId MASTER = 0;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger );  //!< logger for this IO class
};

}  // namespace lama

}  // namespace scai

