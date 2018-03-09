#include <scai/lama/Vector.hpp>
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/solver/Solver.hpp>

/**
 * NewtonStepCG computes search directions for the truncated Newton's method,
 * as used in optimization.
 *
 * TODO: Docs!
 */
template<typename ValueType>
class NewtonStepCG
{
public:

    NewtonStepCG( const scai::lama::Matrix<ValueType>* hessian);

    ValueType getForcingTerm() const;
    void setForcingTerm(const ValueType forcingTerm);

    scai::IndexType getMaxIterations() const;
    void setMaxIterations(scai::IndexType maxIter);

    scai::solver::SolverPtr<ValueType> getPreconditioner() const;
    void setPreconditioner(scai::solver::SolverPtr<ValueType> hessianPreconditioner);

    /**
     * Computes an appropriate search direction for the current gradient and Hessian.
     */
    scai::IndexType computeStep(scai::lama::Vector<ValueType> & dx, const scai::lama::Vector<ValueType> & gradient) const;

private:

    scai::lama::Vector<ValueType> & precondition(const scai::lama::Vector<ValueType> & r) const;

    const scai::lama::Matrix<ValueType> * mHessian;
    ValueType mForcingTerm;
    scai::IndexType mMaxIterations;
    scai::solver::SolverPtr<ValueType> mHessianPreconditioner;

    // Workspace variables, must be pointers as Vector<ValueType> is abstract

    std::unique_ptr<scai::lama::Vector<ValueType>> rPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> zPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> pPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> HpPtr;
    std::unique_ptr<scai::lama::Vector<ValueType>> rMinusGradientPtr;
};
