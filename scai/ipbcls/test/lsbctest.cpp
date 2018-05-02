
#include <scai/ipbcls/ConstrainedLeastSquares.hpp>
#include <scai/ipbcls/test/ConstrainedLeastSquaresProblem.hpp>

#include <cstdlib>

using scai::IndexType;
using scai::lama::DenseVector;
using scai::common::CompareOp;
using scai::ipbcls::ConstrainedLeastSquares;

template <typename Scalar>
bool inInterior( const DenseVector<Scalar> & x, const DenseVector<Scalar> l, const DenseVector<Scalar> u )
{
    return x.all( CompareOp::GT, l ) && x.all( CompareOp::LT, u );
}

typedef scai::DefaultReal ValueType;

int main()
{
    // TODO: Take seed as an optional command-line parameter for
    // reproducible testing.
    using namespace scai::lama;

    const int maxDimension = 100;
    const int maxAttempts = 200;
    const ValueType objTolerance = 1e-4;
    const ValueType resTolerance = 1e-6;

    for ( int i = 0; i < maxAttempts; ++i )
    {
        const IndexType m = rand() % maxDimension;
        const IndexType n = rand() % maxDimension;
        ConstrainedLeastSquaresProblem<ValueType> problem = generateConstrainedLeastSquaresProblem<ValueType>( m, n );

        std::cout << "Solving problem with dimensions " << m << " x " << n << std::endl;

        auto solution = fill<DenseVector<ValueType>>( n, 0.0 );
        ConstrainedLeastSquares<ValueType> solver( problem.A );
        solver.setObjectiveTolerance( objTolerance );
        // solver.setInnerSolverType( InnerSolverType::StandardCG );

        // TODO: Test with non-zero matrix lower bound?
        solver.setResidualTolerance( resTolerance );

        solver.solve( solution, problem.b, problem.l, problem.u );

        if ( !inInterior( solution, problem.l, problem.u ) )
        {
            std::cerr << "Solution x is not in the interior!";
            return 1;
        }

        const auto r_exact  = eval<DenseVector<ValueType>>( problem.A * problem.x - problem.b );
        const auto r_approx = eval<DenseVector<ValueType>>( problem.A * solution - problem.b );

        const auto r_exact_norm = r_exact.l2Norm();
        const auto r_approx_norm = r_approx.l2Norm();

        if ( r_approx_norm > r_exact_norm + objTolerance * r_exact_norm
                && r_approx_norm > resTolerance * problem.b.l2Norm() )
        {
            std::cerr << "Solution does not satisfy accuracy requirements." << std::endl;
            std::cerr << "Exact: " << r_exact_norm << " . Approx: " << r_approx_norm << std::endl;
            return 1;
        }
    }

    return 0;
}
