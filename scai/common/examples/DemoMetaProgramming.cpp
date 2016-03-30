/**
 * @file common/examples/DemoMetaProgramming.cpp
 *
 * @brief Examples for using Template Type lists
 *
 * @author Thomas Brandes
 * @date 15.03.2016
 */

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

#include <iostream>
#include <typeinfo>

// define the list of types

#define MY_TYPELIST TYPELIST( ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )

using namespace scai;
using namespace common;

// This method is instantiated and called for all types of MY_TYPELIST

template<typename T>
void output()
{
    bool isCmp = scalar::isComplex( TypeTraits<T>::stype );

    std::cout << "Running output<" << TypeTraits<T>::id() << ">, is complex = " << isCmp << std::endl;
    T x = 1;
    x = x / 3;
    std::cout << TypeTraits<T>::stype << " : " << x << std::endl;
}

// general defintion of Calling for template list

template<typename TList> struct Calling;

// termination call

template<> struct Calling<mepr::NullType>
{
    static void call()
    {
    }
};

// call output for header T and recursive call for tail of list

template<typename HeadType, typename TailTypes>
struct Calling<mepr::TypeList<HeadType, TailTypes> >
{
    static void call()
    {
        output<HeadType>();
        Calling<TailTypes>::call();
    }
};

int main()
{
    Calling<MY_TYPELIST>::call();
}
