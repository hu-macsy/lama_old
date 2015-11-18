/**
 * @file DynRoutine.hpp
 *
 * @brief Definition of a base class inclusive factory for example of using a library module.
 *
 * @author Thomas Brandes
 * @date 04.11.2015
 */

#include <scai/common/Factory.hpp>

using namespace scai::common;

/** Base class that provides by deriving from Factory a factory with a create routine.
 *
 *  The input value for the create routine is a string, the output value a pointer to
 *  a new created DynRoutine object.
 *
 *  \code
 *       *DynRoutine create( std::string )
 *  \endcode
 */

class DynRoutine  : public scai::common::Factory<std::string, DynRoutine*> 
{
public:

    /** This routine must be provided by all derived classes. */

    virtual void doIt() = 0;

    virtual ~DynRoutine()
    {
    }
};
