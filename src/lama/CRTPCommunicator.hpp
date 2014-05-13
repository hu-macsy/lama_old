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
    
    */

    /** Implementation of pure routine Communicator::swap<double> */

    virtual void swap( double val[], const IndexType n, PartitionId partner ) const
    {
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );
    }

    /** Implementation of pure routine Communicator::swap<float> */

    virtual void swap( float val[], const IndexType n, PartitionId partner ) const
    {
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );
    }

    /** Implementation of pure routine Communicator::swap<IndexType> */

    virtual void swap( IndexType val[], const IndexType n, PartitionId partner ) const
    {
        static_cast<const Derived*>( this )->swapImpl( val, n, partner );
    }

    /***************************************************************
     *  min                                                        *
     **************************************************************/

    virtual IndexType min( const IndexType value ) const
    {
        return static_cast<const Derived*>( this )->minImpl( value );
    }
    
    virtual float min( const float value ) const
    {
        return static_cast<const Derived*>( this )->minImpl( value );
    }
    
    virtual double min( const double value ) const
    {
        return static_cast<const Derived*>( this )->minImpl( value );
    }

    /***************************************************************
     *  max                                                        *
     **************************************************************/

    virtual IndexType max( const IndexType value ) const
    {
        return static_cast<const Derived*>( this )->maxImpl( value );
    }

    virtual double max( const double value ) const
    {
        return static_cast<const Derived*>( this )->maxImpl( value );
    }
    
    virtual float max( const float value ) const
    {
        return static_cast<const Derived*>( this )->maxImpl( value );
    }
    
    /***************************************************************
     *  sum                                                        *
     **************************************************************/

    virtual float sum( const float value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

    virtual double sum( const double value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

    virtual IndexType sum( const IndexType value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

    virtual size_t sum( const size_t value ) const
    {
        return static_cast<const Derived*>( this )->sumImpl( value );
    }

    /***************************************************************
     *  bcast                                                      *
     **************************************************************/

    /** Broadcast an array of doubles from root to all processors in simulation */

    virtual void bcast( double val[], const IndexType n, const PartitionId root ) const
    {
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);
    }

    /** Broadcast an array of IndexType from root to all processors in simulation */

    virtual void bcast( IndexType val[], const IndexType n, const PartitionId root ) const
    {
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);
    }

    /** Broadcast an array of float from root to all processors in simulation */

    virtual void bcast( float val[], const IndexType n, const PartitionId root ) const
    {
        static_cast<const Derived*>( this )->bcastImpl( val, n, root);
    }

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

    /***************************************************************
     *  maxloc                                                     *
     **************************************************************/

    virtual void maxloc( IndexType& val, int& location, PartitionId root ) const
    {
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );
    }

    virtual void maxloc( float& val, int& location, PartitionId root ) const
    {
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );
    }

    virtual void maxloc( double& val, int& location, PartitionId root ) const
    {
        static_cast<const Derived*>( this )->maxlocImpl( val, location, root );
    }

    /***************************************************************
     *  gather                                                     *
     **************************************************************/

    /* Derived classes must provide this routine:

       template<typename T>
       void gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const;
    */

    virtual void gather( IndexType allvals[], const IndexType n, const PartitionId root, const IndexType myvals[] ) const
    {
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );
    }

    virtual void gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const
    {
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );
    }

    virtual void gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const
    {
        static_cast<const Derived*>( this )->gatherImpl( allvals, n, root, myvals );
    }

    /***************************************************************
     *  gatherV                                                    *
     **************************************************************/

    /* Derived classes must provide this routine:

       template<typename T>
       void gatherImpl( T allvals[], const IndexType n, const PartitionId root,
                        const T myvals[], const IndexType sizes[] ) const;
    */

    virtual void gatherV( IndexType allvals[], const IndexType n, const PartitionId root, 
                          const IndexType myvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root, myvals, sizes );
    }

    virtual void gatherV( float allvals[], const IndexType n, const PartitionId root, 
                          const float myvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root, myvals, sizes );
    }

    virtual void gatherV( double allvals[], const IndexType n, const PartitionId root, 
                          const double myvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->gatherVImpl( allvals, n, root, myvals, sizes );
    }

    /***************************************************************
     *  scatter                                                    *
     **************************************************************/

    /* Derived classes must provide this routine:

       template<typename T>
       void scatterImpl( T myvals[], const IndexType n, const PartitionId root, 
                         const T allvals[] ) const;
    */

    virtual void scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const
    {
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );
    }

    virtual void scatter( double myvals[], const IndexType n, const PartitionId root, const double allvals[] ) const
    {
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );
    }

    virtual void scatter( IndexType myvals[], const IndexType n, const PartitionId root, const IndexType allvals[] ) const
    {
        static_cast<const Derived*>( this )->scatterImpl( myvals, n, root, allvals );
    }

    /***************************************************************
     *  scatterV                                                   *
     **************************************************************/

    /* Derived classes must provide this routine:

       template<typename T>
       void scatterVImpl( T myvals[], const IndexType n, const PartitionId root, 
                          const T allvals[], const IndexType sizes[] ) const;
    */

    virtual void scatterV( IndexType myvals[], const IndexType n, const PartitionId root,
                           const IndexType allvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root, allvals, sizes );
    }

    virtual void scatterV( float myvals[], const IndexType n, const PartitionId root,
                           const float allvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root, allvals, sizes );
    }

    virtual void scatterV( double myvals[], const IndexType n, const PartitionId root,
                           const double allvals[], const IndexType sizes[] ) const
    {
        static_cast<const Derived*>( this )->scatterVImpl( myvals, n, root, allvals, sizes );
    }

    /***************************************************************
     *  shift                                                      *
     **************************************************************/

    virtual IndexType shiftData( IndexType recvVals[], const IndexType recvSize,
                                 const IndexType sendVals[], const IndexType sendSize, const int direction ) const
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
 
    virtual IndexType shiftData(
        float recvVals[],
        const IndexType recvSize,
        const float sendVals[],
        const IndexType sendSize,
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
 
    virtual IndexType shiftData(
        double recvVals[],
        const IndexType recvSize,
        const double sendVals[],
        const IndexType sendSize,
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

    /***************************************************************
     *  shiftAsync                                                 *
     **************************************************************/

    virtual SyncToken* shiftAsyncData(
        IndexType recvVals[],
        const IndexType sendVals[],
        const IndexType size,
        const int direction ) const
    {
        if ( direction % getSize() == 0 )
        {
            return defaultShiftAsync( recvVals, sendVals, size, 0 );
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        LAMA_LOG_DEBUG( Derived::logger, "shiftAsyncData<IndexType>, dest = " << dest 
                                         << ", source = " << source << ", size = " << size )
 
        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

    virtual SyncToken* shiftAsyncData(
        float recvVals[],
        const float sendVals[],
        const IndexType size,
        const int direction ) const
    {
        if ( direction % getSize() == 0 )
        {
            return defaultShiftAsync( recvVals, sendVals, size, 0 );
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        LAMA_LOG_DEBUG( Derived::logger, "shiftAsyncData<float>, dest = " << dest 
                                         << ", source = " << source << ", size = " << size )
 
        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

    virtual SyncToken* shiftAsyncData(
        double recvVals[],
        const double sendVals[],
        const IndexType size,
        const int direction ) const
    {
        if ( direction % getSize() == 0 )
        {
            return defaultShiftAsync( recvVals, sendVals, size, 0 );
        }

        PartitionId dest = getNeighbor( direction );
        PartitionId source = getNeighbor( -direction );

        LAMA_LOG_DEBUG( Derived::logger, "shiftAsyncData<double>, dest = " << dest 
                                         << ", source = " << source << ", size = " << size )
 
        return static_cast<const Derived*>( this )->shiftAsyncImpl( recvVals, source, sendVals, dest, size );
    }

    /***************************************************************
     *  exchangeByPlan                                             *
     **************************************************************/

    virtual void exchangeByPlan( IndexType recvVals[], const CommunicationPlan& recvPlan,
                                 const IndexType sendVals[], const CommunicationPlan& sendPlan ) const
    {
        static_cast<const Derived*>( this )->exchangeByPlanImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

    virtual void exchangeByPlan( float recvVals[], const CommunicationPlan& recvPlan,
                                 const float sendVals[], const CommunicationPlan& sendPlan ) const
    {
        static_cast<const Derived*>( this )->exchangeByPlanImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

    virtual void exchangeByPlan( double recvVals[], const CommunicationPlan& recvPlan,
                                 const double sendVals[], const CommunicationPlan& sendPlan ) const
    {
        static_cast<const Derived*>( this )->exchangeByPlanImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

    /***************************************************************
     *  exchangeByPlanAsync                                        *
     **************************************************************/

    virtual SyncToken* exchangeByPlanAsync(
        IndexType recvVals[],
        const CommunicationPlan& recvPlan,
        const IndexType sendVals[],
        const CommunicationPlan& sendPlan ) const
    {
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

    virtual SyncToken* exchangeByPlanAsync(
        float recvVals[],
        const CommunicationPlan& recvPlan,
        const float sendVals[],
        const CommunicationPlan& sendPlan ) const
    {
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

    virtual SyncToken* exchangeByPlanAsync(
        double recvVals[],
        const CommunicationPlan& recvPlan,
        const double sendVals[],
        const CommunicationPlan& sendPlan ) const
    {
        return static_cast<const Derived*>( this )->exchangeByPlanAsyncImpl( recvVals, recvPlan, sendVals, sendPlan );
    }

protected:

    // Default constructor can only be called by derived classes.

    CRTPCommunicator<Derived>( const std::string& type ) : Communicator( type )
    {
    }

};

} // namespace lama

#endif // LAMA_CRTP_MATRIX_STORAGE_HPP_
