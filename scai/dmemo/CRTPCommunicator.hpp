/**
 * @file CRTPCommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief CRTP class that provides the polymorphism for the virtual routines and
 *        so derived classes have to provide only template routines
 * @author Thomas Brandes
 * @date 12.05.2014
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Communicator.hpp>

// internal scai libraris
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/loop.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

/** This template class supports static polymorphism to define
 *  common routines for derived classes of Communicator
 *
 */

template<class Derived>
class COMMON_DLL_IMPORTEXPORT CRTPCommunicator: public Communicator
{
public:

    /* Derived classes must provide this template routine:

     template<typename T>
     void swapImpl( T val[], const IndexType n, PartitionId partner ) const;

     template<typename T>
     T minImpl( const T value ) const

     template<typename T>
     T maxImpl( const T value ) const

     template<typename T>
     T sumImpl( const T value ) const

     template<typename T>
     void bcastImpl( T val[], const IndexType n, const PartitionId root ) const

     template<typename T>
     void gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const;

     template<typename T>
     void gatherVImpl( T allvals[], const IndexType n, const PartitionId root,
     const T myvals[], const IndexType sizes[] ) const;

     template<typename T>
     void scatterImpl( T myvals[], const IndexType n, const PartitionId root,
     const T allvals[] ) const;

     template<typename T>
     void scatterVImpl( T myvals[], const IndexType n, const PartitionId root,
     const T allvals[], const IndexType sizes[] ) const;
     */

#define SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS( _type )                                       \
                                                                                            \
    virtual _type min( const _type value ) const                                            \
    {                                                                                       \
        return static_cast<const Derived*>( this )->minImpl( value );                       \
    }                                                                                       \
                                                                                            \
    virtual void swap(                                                                      \
            _type val[],                                                                    \
            const IndexType n,                                                              \
            PartitionId partner ) const                                                     \
    {                                                                                       \
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );                   \
    }                                                                                       \
                                                                                            \
    virtual _type max( const _type value ) const                                            \
    {                                                                                       \
        return static_cast<const Derived*>( this )->maxImpl( value );                       \
    }                                                                                       \
                                                                                            \
    virtual _type sum( const _type value ) const                                            \
    {                                                                                       \
        return static_cast<const Derived*>( this )->sumImpl( value );                       \
    }                                                                                       \
                                                                                            \
    virtual void bcast( _type val[], const IndexType n, const PartitionId root ) const      \
    {                                                                                       \
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);                      \
    }                                                                                       \
                                                                                            \
    virtual void all2allv(_type* recvVal[],IndexType recvCount[],                           \
                          _type* sendVal[],                                                 \
                          IndexType sendCount[]) const                                      \
    {                                                                                       \
        return static_cast<const Derived*>( this )->all2allvImpl(recvVal,recvCount,         \
                sendVal,sendCount);                                                         \
    }                                                                                       \
                                                                                            \
    virtual void maxloc( _type& val, IndexType& location, PartitionId root ) const          \
    {                                                                                       \
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );             \
    }                                                                                       \
                                                                                            \
    /***************************************************************                        \
         *  gather                                                     *                    \
         **************************************************************/                    \
                                                                                            \
    virtual void gather( _type allvals[], const IndexType n,                                \
                         const PartitionId root, const _type myvals[] ) const               \
    {                                                                                       \
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );        \
    }                                                                                       \
                                                                                            \
    virtual void gatherV(                                                                   \
            _type allvals[], const IndexType n, const PartitionId root,                     \
            const _type myvals[], const IndexType sizes[] ) const                           \
    {                                                                                       \
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root,                 \
                myvals, sizes );                                                            \
    }                                                                                       \
                                                                                            \
    virtual void scatter( _type myvals[],                                                   \
                          const IndexType n,                                                \
                          const PartitionId root,                                           \
                          const _type allvals[] ) const                                     \
    {                                                                                       \
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );       \
    }                                                                                       \
                                                                                            \
    virtual void scatterV(                                                                  \
            _type myvals[],                                                                 \
            const IndexType n,                                                              \
            const PartitionId root,                                                         \
            const _type allvals[],                                                          \
            const IndexType sizes[] ) const                                                 \
    {                                                                                       \
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root,                 \
                allvals, sizes );                                                           \
    }                                                                                       \
                                                                                            \
    virtual IndexType shiftData(                                                            \
            _type recvVals[],                                                               \
            const IndexType recvSize,                                                       \
            const _type sendVals[],                                                         \
            const IndexType sendSize,                                                       \
            const int direction ) const                                                     \
    {                                                                                       \
        return this->shiftDataT( recvVals, recvSize, sendVals, sendSize, direction );       \
    }                                                                                       \
                                                                                            \
    virtual tasking::SyncToken* shiftDataAsync(                                             \
            _type recvVals[],                                                               \
            const _type sendVals[],                                                         \
            const IndexType size,                                                           \
            const int direction ) const                                                     \
    {                                                                                       \
        return this->shiftDataAsyncT( recvVals, sendVals, size, direction );                \
    }                                                                                       \
    virtual void exchangeByPlan( _type recvVals[],                                          \
                                 const CommunicationPlan& recvPlan,                         \
                                 const _type sendVals[],                                    \
                                 const CommunicationPlan& sendPlan ) const                  \
    {                                                                                       \
        static_cast<const Derived*>( this )->exchangeByPlanImpl( recvVals, recvPlan,        \
                sendVals, sendPlan );                                                       \
    }                                                                                       \
                                                                                            \
    /***************************************************************                        \
         *  exchangeByPlanAsync                                        *                    \
         **************************************************************/                    \
                                                                                            \
    virtual tasking::SyncToken* exchangeByPlanAsync(                                        \
            _type recvVals[],                                                               \
            const CommunicationPlan& recvPlan,                                              \
            const _type sendVals[],                                                         \
            const CommunicationPlan& sendPlan ) const                                       \
    {                                                                                       \
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl(                \
                recvVals, recvPlan, sendVals, sendPlan );                                   \
    }

    // define communicator methods for all supported data types

    SCAI_COMMON_LOOP( SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS, SCAI_ALL_TYPES )

#undef SCAI_DMEMO_CRTP_COMMUNICATOR_METHODS

    virtual void bcast( char val[], const IndexType n, const PartitionId root ) const
    {
        static_cast<const Derived*>( this )->bcastImpl( val, n, root );
    }

    /** Additional methods */

    virtual size_t sum( const size_t value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

protected:

    // Default constructor can only be called by derived classes.

    CRTPCommunicator<Derived>( const CommunicatorKind& type )
        : Communicator( type )
    {
    }

private:

    template<typename T>
    IndexType shiftDataT(
        T recvVals[],
        const IndexType recvSize,
        const T sendVals[],
        const IndexType sendSize,
        const int direction ) const
    {
        SCAI_LOG_DEBUG( Derived::logger,
                        *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

        if ( direction % getSize() == 0 )
        {
            return shift0( recvVals, recvSize, sendVals, sendSize );
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );
        return static_cast<const Derived*>( this )->shiftImpl( recvVals, recvSize, source, sendVals, sendSize, dest );
    }

    template<typename T>
    tasking::SyncToken* shiftDataAsyncT( T recvVals[], const T sendVals[], const IndexType size, const int direction ) const
    {
        if ( direction % getSize() == 0 )
        {
            shift0( recvVals, size, sendVals, size );
            return new tasking::NoSyncToken();
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        SCAI_LOG_DEBUG( Derived::logger,
                        "shiftDataAsync<" << typeid( T ).name() << ">, dest = " << dest << ", source = " << source << ", size = " << size )

        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

};

} /* end namespace dmemo */

} /* end namespace scai */
