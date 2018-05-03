#include "NewtonStepCG.hpp"

#include <scai/tracing.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/dmemo/Distribution.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/instantiate.hpp>

using scai::IndexType;
using scai::lama::Matrix;
using scai::lama::DenseVector;
using scai::lama::Vector;

namespace
{
}

template<typename ValueType>
NewtonStepCG<ValueType>::NewtonStepCG( const Matrix<ValueType>* hessian )
    : mHessian( hessian ),
      mForcingTerm( 0.5 ),
      mMaxIterations( 0 ),
      mHessianPreconditioner( new scai::solver::TrivialPreconditioner<ValueType>( "Identitiy preconditioner" ) )
{
    rPtr.reset( hessian->newTargetVector() );
    zPtr.reset( hessian->newSourceVector() );
    pPtr.reset( hessian->newSourceVector() );
    HpPtr.reset( hessian->newTargetVector() );
    rMinusGradientPtr.reset( hessian->newTargetVector() );

    mHessianPreconditioner->initialize(*hessian);
}

template<typename ValueType>
ValueType NewtonStepCG<ValueType>::getForcingTerm() const
{
    return mForcingTerm;
}

template<typename ValueType>
scai::solver::SolverPtr<ValueType> NewtonStepCG<ValueType>::getPreconditioner() const
{
    return mHessianPreconditioner;
}

template<typename ValueType>
IndexType NewtonStepCG<ValueType>::getMaxIterations() const
{
    return mMaxIterations;
}

template<typename ValueType>
void NewtonStepCG<ValueType>::setForcingTerm( ValueType forcingTerm )
{
    SCAI_ASSERT_ERROR( forcingTerm < 1 && forcingTerm > 0, "Forcing term " << forcingTerm << " must be in the open interval (0, 1)" )
    mForcingTerm = forcingTerm;
}

template<typename ValueType>
void NewtonStepCG<ValueType>::setMaxIterations( IndexType maxIter )
{
    SCAI_ASSERT_ERROR( maxIter >= 0, "Maximum number of iterations must be non-negative." )
    mMaxIterations = maxIter;
}

template<typename ValueType>
void NewtonStepCG<ValueType>::setPreconditioner( scai::solver::SolverPtr<ValueType> hessianPreconditioner )
{
    mHessianPreconditioner = hessianPreconditioner;
}

template<typename ValueType>
Vector<ValueType>& NewtonStepCG<ValueType>::precondition( const Vector<ValueType>& r ) const
{
    // This method exists to wrap the zero-initialization together with the preconditioning.
    // Ideally, we wouldn't need to zero-initialize z to store the preconditioning result,
    // but because the preconditioner is actually a solver, and not just a linear operator,
    // it assumes that z is an "initial guess" (which it is not in our case). Hence we're
    // better off zero-initializing it first.

    *zPtr = 0;
    getPreconditioner()->solve( *zPtr, r );
    return *zPtr;
}

template<typename ValueType>
IndexType NewtonStepCG<ValueType>::computeStep( Vector<ValueType>& dx, const Vector<ValueType>& gradient ) const
{
    SCAI_REGION( "NewtonStepCG::computeStep" );

    const ValueType theta = getForcingTerm();

    const Matrix<ValueType> & H = *mHessian;

    Vector<ValueType>& r = *rPtr;
    Vector<ValueType>& z = *zPtr;
    Vector<ValueType>& p = *pPtr;
    Vector<ValueType>& Hp = *HpPtr;
    Vector<ValueType>& rMinusGradient = *rMinusGradientPtr;

    ValueType quadraticObjective = 1;
    ValueType rTz;

    // Use initial guess dx only if it is a descent direction for the
    // quadratic model. Otherwise, set dx = 0.
    do
    {
        // r = -g - H dx
        r = -1 * gradient;
        r -= H * dx;
        z = precondition( r );
        p = z;
        Hp = H * p;
        rTz = r.dotProduct( z );

        if ( rTz == 0 )
        {
            SCAI_ASSERT_EQ_ERROR( r.l2Norm(), ValueType( 0 ), "Preconditioner is not positive definite." );

            // If r == 0, then we're already arrived at the solution. We break here because
            // we'll otherwise end up dividing by zero later (since p, Hp are both zero as well).
            return IndexType( 0 );
        }

        rMinusGradient = r - gradient;
        quadraticObjective = -0.5 * dx.dotProduct( rMinusGradient );

        if ( quadraticObjective > 0 )
        {
            dx = 0;
        }
    }
    while ( quadraticObjective > 0 );

    IndexType k = 1;

    while ( k <= getMaxIterations() || getMaxIterations() == IndexType( 0 ) )
    {
        const ValueType alpha = rTz / p.dotProduct( Hp );
        dx = dx + alpha * p;
        r = r - alpha * Hp;

        rMinusGradient = r - gradient;
        const ValueType newQuadraticObjective = - 0.5 * dx.dotProduct( rMinusGradient );

        if ( newQuadraticObjective - quadraticObjective >= theta * newQuadraticObjective / ValueType ( k ) )
        {
            // See the very short paper
            // SG Nash, A Sofer: "Assessing a search direction within a truncated-Newton method"
            // for a justification for this stopping criterion. Essentially, it is a far better predictor
            // for the quality of a search direction in truncated Newton methods than
            // residual-based criteria.
            break;
        }

        z = precondition( r );
        const ValueType rTzNew = r.dotProduct( z );

        if ( rTzNew == ValueType( 0 ) ) 
        {
            // might happen if solution is found directly.
            break;
        }
        const ValueType beta = rTzNew / rTz;
        p = z + beta * p;

        quadraticObjective = newQuadraticObjective;
        rTz = rTzNew;
        Hp = H * p;
        k += 1;
    }

    return k - 1;
}

// Template instantiation only for real (non-complex) types

SCAI_COMMON_INST_CLASS( NewtonStepCG, SCAI_REAL_TYPES_HOST )
