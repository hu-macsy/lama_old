/*
 * XDRUtils.hpp
 *
 *  Created on: Mar 8, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/mepr/TypeList.hpp>

namespace scai {

namespace lama {

namespace mepr {

/*
 * _RegistratorVO
 */

template<typename ValueType, typename TList> struct IOWrapper;

template<typename ValueType>
struct IOWrapper<ValueType, common::mepr::NullType>
{
    static void readXDR( const common::scalar::ScalarType, XDRFileStream&, ValueType*, const IndexType ) {}
    static void writeXDR( const common::scalar::ScalarType, XDRFileStream&, const ValueType*, const IndexType ) {}
    static void readBinary( const common::scalar::ScalarType, std::fstream&, ValueType*, const IndexType ) {}
    static void writeBinary( const common::scalar::ScalarType, std::fstream&, const ValueType*, const IndexType ) {}
};

template<typename ValueType, typename H, typename T>
struct IOWrapper<ValueType, common::mepr::TypeList<H,T> >
{
    static void readXDR( const common::scalar::ScalarType dataType, XDRFileStream& in, ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::readXDR<H,ValueType>( in, data, n );
        }
        else
        {
            IOWrapper<ValueType, T>::readXDR( dataType, in, data, n );
        }
    }

    static void writeXDR( const common::scalar::ScalarType dataType, XDRFileStream& out, const ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::writeXDR<H,ValueType>( out, data, n );
        }
        else
        {
            IOWrapper<ValueType, T>::writeXDR( dataType, out, data, n );
        }
    }

    static void readBinary( const common::scalar::ScalarType dataType, std::fstream& in, ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            FileIO::readBinaryData<H, ValueType>( in, data, n );
        }
        else
        {
            IOWrapper<ValueType, T>::writeBinary( dataType, in, data, n );
        }
    }

    static void writeBinary( const common::scalar::ScalarType dataType, std::fstream& out, const ValueType* data, const IndexType n )
    {
        if( dataType == common::getScalarType<H>() )
        {
            IOUtils::writeBinary<H,ValueType>( out, data, n );
        }
        else
        {
            IOWrapper<ValueType, T>::writeBinary( dataType, out, data, n );
        }
    }
};

/*
 * RegistratorVO
 */

//template<template<typename,typename> class R, typename ValueType, typename TList1, typename TList2> struct XDRWrapper;
//
//template<template<typename,typename> class R, typename ValueType >
//struct XDRWrapper<R, ValueType, common::mepr::NullType, common::mepr::NullType >
//{
//    static void read( const common::scalar::ScalarType, XDRFileStream& , ValueType[], const IndexType ) {}
//};
//
//template<template<typename,typename> class R, typename ValueType, typename TList >
//struct XDRWrapper<R, ValueType, common::mepr::NullType, TList >
//{
//    static void read( const common::scalar::ScalarType, XDRFileStream& , ValueType[], const IndexType ) {}
//};
//
//template<template<typename,typename> class R, typename ValueType, typename H1, typename T1, typename H2, typename T2>
//struct XDRWrapper< R, ValueType, common::mepr::TypeList<H1,T1>, common::mepr::TypeList<H2, T2> >
//{
//    static void read( const common::scalar::ScalarType dataType, XDRFileStream& in, ValueType data[], const IndexType n )
//    {
//        IOWrapper<R, ValueType, H1, common::mepr::TypeList<H2,T2> >::read( dataType, in, data, n );
//
//        XDRWrapper<R, ValueType, T1, common::mepr::TypeList<H2,T2> >::read( dataType, in, data, n );
//    }
//};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
