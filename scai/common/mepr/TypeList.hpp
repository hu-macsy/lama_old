#pragma once

namespace scai {

namespace common {

namespace mepr {

class NullType {};

template<typename H, typename T>
struct TypeList
{
    typedef H head;
    typedef T tail;
};

#define TYPELIST_1( T1 ) common::mepr::TypeList<T1,common::mepr::NullType>
#define TYPELIST_2( T1, T2 ) common::mepr::TypeList<T1,TYPELIST_1( T2 ) >
#define TYPELIST_3( T1, T2, T3 ) common::mepr::TypeList<T1,TYPELIST_2( T2, T3 ) >
#define TYPELIST_4( T1, T2, T3, T4 ) common::mepr::TypeList<T1,TYPELIST_3( T2, T3, T4 ) >
#define TYPELIST_5( T1, T2, T3, T4, T5 ) common::mepr::TypeList<T1,TYPELIST_4( T2, T3, T4, T5 ) >
#define TYPELIST_6( T1, T2, T3, T4, T5, T6 ) common::mepr::TypeList<T1,TYPELIST_5( T2, T3, T4, T5, T6 ) >
#define TYPELIST_7( T1, T2, T3, T4, T5, T6, T7 ) common::mepr::TypeList<T1,TYPELIST_6( T2, T3, T4, T5, T6, T7 ) >
#define TYPELIST_8( T1, T2, T3, T4, T5, T6, T7, T8 ) common::mepr::TypeList<T1,TYPELIST_7( T2, T3, T4, T5, T6, T7, T8 ) >

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
