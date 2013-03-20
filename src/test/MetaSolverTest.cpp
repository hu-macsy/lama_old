/**
 * @file MetaSolverTest.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief Contains the implementation of the class JacobiTest.cpp
 * @author: Kai Buschulte
 * @date 07.05.2012
 * $
 **/

#include <test/MetaSolverTest.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/SimpleAMG.hpp>
#include <lama/solver/SOR.hpp>
#include <lama/solver/GMRES.hpp>
#include <lama/solver/InverseSolver.hpp>
#include <lama/solver/DefaultJacobi.hpp>
#include <lama/solver/SpecializedJacobi.hpp>

#include <lama/solver/logger/OpenMPTimer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/Timer.hpp>

#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/criteria/ResidualStagnation.hpp>

#include <lama/solver/MetaSolver.hpp>
#include <lama/solver/creator/CriteriaCreator.hpp>

#include <lama/matutils/MatrixCreator.hpp>
#include <lama/DenseVector.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <lama/norm/L2Norm.hpp>

#include <test/TestMacros.hpp>

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MetaSolverTest )
;

using namespace boost;
using namespace lama;

LAMA_LOG_DEF_LOGGER( logger, "Test.MetaSolverTest" );

/* --------------------------------------------------------------------- */

MetaSolverTestImpl::MetaSolverTestImpl( std::string& config, SolverPtr native )
    : mConfiguration( config ), mNative( native )
{
    const IndexType N1 = 4;
    const IndexType N2 = 4;

    LAMA_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );

    LAMA_LOG_INFO( logger, "Defining Matrix, solution Vector and rhs." );

    CSRSparseMatrix<ValueType>* csrMatrix = new CSRSparseMatrix<ValueType>();

    MatrixCreator<ValueType>::buildPoisson2D( *csrMatrix, 9, N1, N2 );
    mMatrix.reset( csrMatrix );

    mSolution.reset( new DenseVector<ValueType>( mMatrix->getColDistributionPtr(), 1.5 ) );
    mRhs.reset( new DenseVector<ValueType>( mMatrix->getDistributionPtr(), 0.5 ) );
}

MetaSolverTestImpl::~MetaSolverTestImpl()
{
}

/* --------------------------------------------------------------------- */

void MetaSolverTestImpl::configTest()
{
    MetaSolver* metaSolver = new MetaSolver( "MetaSolverTest" );

    LAMA_LOG_INFO( logger, "Created MetaSolver instance. Parsing '" << mConfiguration << "' now." );

    metaSolver->parseConfiguration( mConfiguration );

    LAMA_LOG_INFO( logger, "Parsed Configuration." );

    mMetaSolver.reset( metaSolver );
}

/* --------------------------------------------------------------------- */

void MetaSolverTestImpl::initializeTest()
{
    LAMA_LOG_INFO( logger, "Initializing MetaSolver with Matrix " << *mMatrix << "." );
    mMetaSolver->initialize( *mMatrix );
}

/* --------------------------------------------------------------------- */

void MetaSolverTestImpl::solveTest()
{
    LAMA_LOG_INFO( logger, "MetaSolver solve with solution vector " << *mSolution << " and rhs " << *mRhs );

    mMetaSolver->solve( *mSolution, *mRhs );
}

/* --------------------------------------------------------------------- */

void MetaSolverTestImpl::testCompareNative()
{
    LAMA_LOG_INFO( logger, "Native implementation: Calling initialize and solve." );

    mNative->initialize( *mMatrix );

    DenseVector<ValueType> solution( mSolution->getDistributionPtr(), 1.5 );
    mNative->solve( solution, *mRhs );

    LAMA_LOG_INFO( logger, "Comparing solution of native and MetaSolver implementations." );

    Scalar nativeNorm = l2Norm( solution );
    Scalar metaNorm = l2Norm( *mSolution );
    LAMA_LOG_INFO( logger, "Native norm: " << nativeNorm << ". MetaSolver Norm: " << metaNorm << "." );

    LAMA_CHECK_SCALAR_CLOSE( nativeNorm, metaNorm, double, 1E-3 );
}

/* --------------------------------------------------------------------- */

void MetaSolverTestImpl::runTests()
{
    configTest();
    initializeTest();
    solveTest();
    testCompareNative();
}

/* ---------------------Test Scenarios-defined here-------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithSOR )
{
    std::string configuration = "SOR root{omega=1.5;}";
    SOR* sor = new SOR( "SOR_Native" );
    sor->setOmega( Scalar( 1.5 ) );
    SolverPtr nativeImpl( sor );
    MetaSolverTestImpl metaSolverTestWithSOR( configuration, nativeImpl );
    metaSolverTestWithSOR.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithCG )
{
    std::string configuration = "CG cg{} CG root{preconditioner=cg;}";
    SolverPtr cg( new CG( "CG_Native_Preconditioner" ) );
    CG* root = new CG( "CG_Native" );
    root->setPreconditioner( cg );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithCG( configuration, nativeImpl );
    metaSolverTestWithCG.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithGMRES )
{
    std::string configuration = "GMRES root{krylovDim=3;}";
    GMRES* root = new GMRES( "GMRES_Native" );
    root->setKrylovDim( 3 );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithGMRES( configuration, nativeImpl );
    metaSolverTestWithGMRES.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithInversSolver )
{
    std::string configuration = "InverseSolver root{}";
    InverseSolver* inv = new InverseSolver( "InverseSolver_Native" );
    SolverPtr nativeImpl( inv );
    MetaSolverTestImpl metaSolverTestWithInv( configuration, nativeImpl );
    metaSolverTestWithInv.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithSimpleAMG )
{
    std::string configuration =
        "SOR sor{} InverseSolver inv{} \
            SimpleAMG root{smoother=sor; coarseLevelSolver=inv; \
            maxLevels=10;minVarsCoarseLevel=100;}";
    SolverPtr sor( new SOR( "Smoother_Native" ) );
    SolverPtr inv( new InverseSolver( "CoarseLevel_Native" ) );
    SimpleAMG* root = new SimpleAMG( "SimpleAMG_Native" );
    root->setSmoother( sor );
    root->setCoarseLevelSolver( inv );
    root->setMaxLevels( 10 );
    root->setMinVarsCoarseLevel( 100 );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithSAMG( configuration, nativeImpl );
    metaSolverTestWithSAMG.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithDefaultJacobi )
{
    std::string configuration = "DefaultJacobi root{omega=1.5;}";
    DefaultJacobi* root = new DefaultJacobi( "DefaultJacobi_Native" );
    root->setOmega( Scalar( 1.5 ) );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithDJ( configuration, nativeImpl );
    metaSolverTestWithDJ.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithSpecializedJacobi )
{
    std::string configuration = "SpecializedJacobi root{omega=1.5;}";
    SpecializedJacobi* root = new SpecializedJacobi( "SpecializedJacobi_Native" );
    root->setOmega( Scalar( 1.5 ) );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithSJ( configuration, nativeImpl );
    metaSolverTestWithSJ.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithStoppingCriteria )
{
    std::string configuration =
        "Criterion testCriteria = ( ( IterationCount( 10 ) OR \
            ResidualThreshold( L2Norm, 0.01, Absolute ) ) OR \
            ResidualStagnation(MaxNorm, 3, 0.001 )); \
            CG root { stoppingCriterion = testCriteria; }";
    CG* root = new CG( "CG_Native" );

    CriterionPtr iterC( new IterationCount( 10 ) );
    CriterionPtr resThr( new ResidualThreshold( NormPtr( new L2Norm() ), 0.01, ResidualThreshold::Absolute ) );
    CriterionPtr resStag( new ResidualStagnation( NormPtr( new MaxNorm() ), 3, 0.001 ) );
    CriterionPtr itResThr( new Criterion( iterC, resThr, Criterion::OR ) );
    CriterionPtr criteria( new Criterion( itResThr, resStag, Criterion::OR ) );

    root->setStoppingCriterion( criteria );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithSOR( configuration, nativeImpl );
    metaSolverTestWithSOR.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( metaSolverTestWithLogger )
{
    std::string configuration =
        "CommonLogger loggerDef( TestLogger, noLogging, toConsoleOnly, OpenMPTimer ); \
            CG root { logger = loggerDef; }";
    LoggerPtr logger(
        new CommonLogger( "TestLogger_Native", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new OpenMPTimer() ) ) );
    CG* root = new CG( "CG_Native" );
    SolverPtr nativeImpl( root );
    MetaSolverTestImpl metaSolverTestWithSOR( configuration, nativeImpl );
    metaSolverTestWithSOR.runTests();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    MetaSolver metaSolver( "metaSolver", "CG root{}" );
    LAMA_WRITEAT_TEST( metaSolver );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
