/**
 * @file CRTPCommunicator.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @since 1.0.1
 */

#ifndef LAMA_CRTP_COMMUNICATOR_HPP_
#define LAMA_CRTP_COMMUNICATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Communicator.hpp>
#include <vector>

namespace lama
{

/** This template class supports static polymorphism to define
 *  common routines for derived classes of Communicator
 *
 */

template<class Derived>
class LAMA_DLL_IMPORTEXPORT CRTPCommunicator: public Communicator
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
    
    /***************************************************************
     *  bcast                                                      *
     **************************************************************/

    /** Broadcast of a string, done in two steps. */

    virtual void bcast( std::string& val, const PartitionId root ) const
    {
        bool isRoot = getRank() == root;

        int len = 0;

        if ( isRoot )
        {
            len = val.length();
        }

        // step 1: broadcast length of string

        static_cast<const Derived*>( this )->bcastImpl( &len, 1, root );

        std::vector<char> buffer( len+1 );

        char* strptr = buffer.data();

        if ( isRoot )
        {
            strptr = const_cast<char *>( val.c_str() );
        }

        static_cast<const Derived*>( this )->bcastImpl( strptr, len + 1, root );

        if ( !isRoot )
        {
            val = strptr;
        }
    }

#define CRTP_COMMUNICATOR_METHODS( TypeId )                                             \
                                                                                        \
    virtual TypeId min( const TypeId value ) const                                      \
    {                                                                                   \
        return static_cast<const Derived*>( this )->minImpl( value );                   \
    }                                                                                   \
                                                                                        \
    virtual void swap( TypeId val[], const IndexType n, PartitionId partner ) const     \
    {                                                                                   \
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );               \
    }                                                                                   \
                                                                                        \
    virtual TypeId max( const TypeId value ) const                                      \
    {                                                                                   \
        return static_cast<const Derived*>( this )->maxImpl( value );                   \
    }                                                                                   \
                                                                                        \
    virtual TypeId sum( const TypeId value ) const                                      \
    {                                                                                   \
        return static_cast<const Derived*>( this )->sumImpl( value );                   \
    }                                                                                   \
                                                                                        \
    virtual void bcast( TypeId val[], const IndexType n, const PartitionId root ) const \
    {                                                                                   \
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);                  \
    }                                                                                   \
                                                                                        \
    virtual void maxloc( TypeId& val, IndexType& location, PartitionId root ) const     \
    {                                                                                   \
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );         \
    }                                                                                   \
                                                                                        \
    /***************************************************************                    \
     *  gather                                                     *                    \
     **************************************************************/                    \
                                                                                        \
    virtual void gather( TypeId allvals[], const IndexType n,                           \
                         const PartitionId root, const TypeId myvals[] ) const          \
    {                                                                                   \
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );    \
    }                                                                                   \
                                                                                        \
    virtual void gatherV( TypeId allvals[], const IndexType n, const PartitionId root,  \
                          const TypeId myvals[], const IndexType sizes[] ) const        \
    {                                                                                   \
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root,             \
                                                          myvals, sizes );              \
    }                                                                                   \
                                                                                        \
    virtual void scatter( TypeId myvals[], const IndexType n, const PartitionId root,   \
                          const TypeId allvals[] ) const                                \
    {                                                                                   \
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );   \
    }                                                                                   \
                                                                                        \
    virtual void scatterV( TypeId myvals[], const IndexType n, const PartitionId root,  \
                           const TypeId allvals[], const IndexType sizes[] ) const      \
    {                                                                                   \
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root,             \
                                                           allvals, sizes );            \
    }                                                                                   \
                                                                                        \
    virtual IndexType shiftData(                                                        \
        TypeId recvVals[],                                                              \
        const IndexType recvSize,                                                       \
        const TypeId sendVals[],                                                        \
        const IndexType sendSize,                                                       \
        const int direction ) const                                                     \
    {                                                                                   \
        return this->shiftDataT( recvVals, recvSize, sendVals, sendSize, direction );   \
    }                                                                                   \
                                                                                        \
    virtual SyncToken* shiftAsyncData(                                                  \
        TypeId recvVals[],                                                              \
        const TypeId sendVals[],                                                        \
        const IndexType size,                                                           \
        const int direction ) const                                                     \
    {                                                                                   \
        return this->shiftAsyncDataT( recvVals, sendVals, size, direction );            \
    }                                                                                   \
                                                                                        \
    virtual void exchangeByPlan( TypeId recvVals[],                                     \
                                 const CommunicationPlan& recvPlan,                     \
                                 const TypeId sendVals[],                               \
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
    virtual SyncToken* exchangeByPlanAsync(                                             \
        TypeId recvVals[],                                                              \
        const CommunicationPlan& recvPlan,                                              \
        const TypeId sendVals[],                                                        \
        const CommunicationPlan& sendPlan ) const                                       \
    {                                                                                   \
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl(            \
                       recvVals, recvPlan, sendVals, sendPlan );                        \
    }

    CRTP_COMMUNICATOR_METHODS( int )
    CRTP_COMMUNICATOR_METHODS( float )
    CRTP_COMMUNICATOR_METHODS( double )

#ifdef LAMA_USE_COMPLEX_FLOAT                   
    CRTP_COMMUNICATOR_METHODS( ComplexFloat )
#endif

#ifdef LAMA_USE_COMPLEX_DOUBLE
    CRTP_COMMUNICATOR_METHODS( ComplexDouble )
#endif

    /** Additional methods */

    virtual size_t sum( const size_t value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

protected:

    // Default constructor can only be called by derived classes.

    CRTPCommunicator<Derived>( const std::string& type ) : Communicator( type )
    {
    }

private:

    template<typename T>
    IndexType shiftDataT( T recvVals[], const IndexType recvSize,  
                          const T sendVals[], const IndexType sendSize,
                          const int direction ) const   
    { 
        LAMA_LOG_DEBUG( Derived::logger, *this << ": shift, direction = " << direction 
                                          << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

        if ( direction % getSize() == 0 )
        {
            return shift0( recvVals, recvSize, sendVals, sendSize );
        }
    
        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );
        return static_cast<const Derived*>( this )->shiftImpl( recvVals, recvSize, source, sendVals, sendSize, dest );
    }
 
    template<typename T>
    SyncToken* shiftAsyncDataT( T recvVals[], const T sendVals[],
                                const IndexType size, const int direction ) const
    {
        if ( direction % getSize() == 0 )
        {
            return defaultShiftAsync( recvVals, sendVals, size, 0 );
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        LAMA_LOG_DEBUG( Derived::logger, "shiftAsyncData<" << typeid(T).name() << ">, dest = " << dest 
                                         << ", source = " << source << ", size = " << size )
 
        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

    };

} // namespace lama

#endif // LAMA_CRTP_MATRIX_STORAGE_HPP_
