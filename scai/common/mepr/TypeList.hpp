#pragma once

namespace scai
{

namespace common
{

namespace mepr
{

class NullType {};

template<typename H, typename T>
struct TypeList
{
    typedef H head;
    typedef T tail;
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */

#define TYPELIST_1( T1 ) scai::common::mepr::TypeList<T1,common::mepr::NullType>
#define TYPELIST_2( T1, T2 ) scai::common::mepr::TypeList<T1,TYPELIST_1( T2 ) >
#define TYPELIST_3( T1, T2, T3 ) scai::common::mepr::TypeList<T1,TYPELIST_2( T2, T3 ) >
#define TYPELIST_4( T1, T2, T3, T4 ) scai::common::mepr::TypeList<T1,TYPELIST_3( T2, T3, T4 ) >
#define TYPELIST_5( T1, T2, T3, T4, T5 ) scai::common::mepr::TypeList<T1,TYPELIST_4( T2, T3, T4, T5 ) >
#define TYPELIST_6( T1, T2, T3, T4, T5, T6 ) scai::common::mepr::TypeList<T1,TYPELIST_5( T2, T3, T4, T5, T6 ) >
#define TYPELIST_7( T1, T2, T3, T4, T5, T6, T7 ) scai::common::mepr::TypeList<T1,TYPELIST_6( T2, T3, T4, T5, T6, T7 ) >
#define TYPELIST_8( T1, T2, T3, T4, T5, T6, T7, T8 ) scai::common::mepr::TypeList<T1,TYPELIST_7( T2, T3, T4, T5, T6, T7, T8 ) >
