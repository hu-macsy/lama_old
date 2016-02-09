/**
 * @file CRTPCommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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

#include <scai/common/preprocessor.hpp>

// std
#include <vector>

namespace scai
{

namespace lama
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

#define CRTP_COMMUNICATOR_METHODS( z, I, _ )                                            \
    \
    virtual ARRAY_TYPE##I min( const ARRAY_TYPE##I value ) const                        \
    {                                                                                   \
        return static_cast<const Derived*>( this )->minImpl( value );                   \
    }                                                                                   \
    \
    virtual void swap(                                                                  \
            ARRAY_TYPE##I val[],                                                            \
            const IndexType n,                                                              \
            PartitionId partner ) const                                                     \
    {                                                                                   \
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );               \
    }                                                                                   \
    \
    virtual ARRAY_TYPE##I max( const ARRAY_TYPE##I value ) const                        \
    {                                                                                   \
        return static_cast<const Derived*>( this )->maxImpl( value );                   \
    }                                                                                   \
    \
    virtual ARRAY_TYPE##I sum( const ARRAY_TYPE##I value ) const                        \
    {                                                                                   \
        return static_cast<const Derived*>( this )->sumImpl( value );                   \
    }                                                                                   \
    \
    virtual void bcast( ARRAY_TYPE##I val[], const IndexType n, const PartitionId root ) const \
    {                                                                                   \
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);                  \
    }                                                                                   \
    \
    virtual void all2allv(ARRAY_TYPE##I* recvVal[],IndexType recvCount[],               \
                        ARRAY_TYPE##I* sendVal[],                                       \
                            IndexType sendCount[]) const                                \
    {                                                                                   \
        return static_cast<const Derived*>( this )->all2allvImpl(recvVal,recvCount,     \
                                                  sendVal,sendCount);                   \
    }                                                                                   \
    \
    virtual void maxloc( ARRAY_TYPE##I& val, IndexType& location, PartitionId root ) const     \
    {                                                                                   \
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );         \
    }                                                                                   \
    \
    /***************************************************************                    \
     *  gather                                                     *                    \
     **************************************************************/                    \
    \
    virtual void gather( ARRAY_TYPE##I allvals[], const IndexType n,                    \
                         const PartitionId root, const ARRAY_TYPE##I myvals[] ) const   \
    {                                                                                   \
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );    \
    }                                                                                   \
    \
    virtual void gatherV(                                                               \
            ARRAY_TYPE##I allvals[], const IndexType n, const PartitionId root,             \
            const ARRAY_TYPE##I myvals[], const IndexType sizes[] ) const                   \
    {                                                                                   \
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root,             \
                myvals, sizes );              \
    }                                                                                   \
    \
    virtual void scatter( ARRAY_TYPE##I myvals[],                                       \
                          const IndexType n,                                            \
                          const PartitionId root,                                       \
                          const ARRAY_TYPE##I allvals[] ) const                         \
    {                                                                                   \
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );   \
    }                                                                                   \
    \
    virtual void scatterV(                                                              \
            ARRAY_TYPE##I myvals[],                                                         \
            const IndexType n,                                                              \
            const PartitionId root,                                                         \
            const ARRAY_TYPE##I allvals[],                                                  \
            const IndexType sizes[] ) const                                                 \
    {                                                                                   \
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root,             \
                allvals, sizes );            \
    }                                                                                   \
    \
    virtual IndexType shiftData(                                                        \
            ARRAY_TYPE##I recvVals[],                                                       \
            const IndexType recvSize,                                                       \
            const ARRAY_TYPE##I sendVals[],                                                 \
            const IndexType sendSize,                                                       \
            const int direction ) const                                                     \
    {                                                                                   \
        return this->shiftDataT( recvVals, recvSize, sendVals, sendSize, direction );   \
    }                                                                                   \
    \
    virtual tasking::SyncToken* shiftDataAsync(                                                  \
            ARRAY_TYPE##I recvVals[],                                                       \
            const ARRAY_TYPE##I sendVals[],                                                 \
            const IndexType size,                                                           \
            const int direction ) const                                                     \
    {                                                                                   \
        return this->shiftDataAsyncT( recvVals, sendVals, size, direction );            \
    }                                                                                   \
    virtual void exchangeByPlan( ARRAY_TYPE##I recvVals[],                              \
                                 const CommunicationPlan& recvPlan,                     \
                                 const ARRAY_TYPE##I sendVals[],                        \
                                 const CommunicationPlan& sendPlan ) const              \
    {                                                                                   \
        static_cast<const Derived*>( this )->exchangeByPlanImpl( recvVals, recvPlan,    \
                sendVals, sendPlan );  \
    }                                                                                   \
    \
    /***************************************************************                    \
     *  exchangeByPlanAsync                                        *                    \
     **************************************************************/                    \
    \
    virtual tasking::SyncToken* exchangeByPlanAsync(                                             \
            ARRAY_TYPE##I recvVals[],                                                       \
            const CommunicationPlan& recvPlan,                                              \
            const ARRAY_TYPE##I sendVals[],                                                 \
            const CommunicationPlan& sendPlan ) const                                       \
    {                                                                                   \
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl(            \
                recvVals, recvPlan, sendVals, sendPlan );                        \
    }

    // define communicator methods for all supported data types

    BOOST_PP_REPEAT( ARRAY_TYPE_CNT, CRTP_COMMUNICATOR_METHODS, _ )

#undef CRTP_COMMUNICATOR_METHODS

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

    CRTPCommunicator<Derived>( const communicator::CommunicatorKind& type )
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

        if( direction % getSize() == 0 )
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
        if( direction % getSize() == 0 )
        {
            shift0( recvVals, size, sendVals, size );
            return new tasking::NoSyncToken();
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        SCAI_LOG_DEBUG( Derived::logger,
                        "shiftDataAsync<" << typeid(T).name() << ">, dest = " << dest << ", source = " << source << ", size = " << size )

        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

};

} /* end namespace lama */

} /* end namespace scai */
