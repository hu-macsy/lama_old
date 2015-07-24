/**
 * @file Richardson.cpp
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
 * @brief Richardson.cpp
 * @author David Schissler
 * @date 17.04.2015
 * @since 
 */

 // hpp
#include <lama/solver/Richardson.hpp>
// others
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/DenseVector.hpp>


namespace lama{

LAMA_LOG_DEF_LOGGER( Richardson::logger, "Solver.Richardson" )

Richardson::Richardson( const std::string& id )
    : OmegaSolver( id,(Scalar) -1.0){}

Richardson::Richardson( const std::string& id, const Scalar omega )
    : OmegaSolver( id, omega ){}

Richardson::Richardson( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id ,(Scalar) -1.0,logger){}

Richardson::Richardson( const std::string& id, const Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger ){}

Richardson::Richardson( const Richardson& other )
    : OmegaSolver( other ){}



Richardson::RichardsonRuntime::RichardsonRuntime()
    : OmegaSolverRuntime(){}

Richardson::~Richardson(){}

Richardson::RichardsonRuntime::~RichardsonRuntime(){}



void Richardson::initialize( const Matrix& coefficients ){
    LAMA_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    Solver::initialize( coefficients );

    if(mOmega == -1.0){
        L2Norm n;
        Scalar bound = 2.0/n.apply( coefficients );
        mOmega = (2.0/3.0 * bound);
    }
}

void Richardson::solveInit( Vector& solution, const Vector& rhs ){
    RichardsonRuntime& runtime = getRuntime();
    //Check if oldSolution already exists, if not create copy of solution
    if ( !runtime.mOldSolution.get() )
        runtime.mOldSolution.reset( solution.clone() );
    
    runtime.mProxyOldSolution = runtime.mOldSolution.get();

    IterativeSolver::solveInit( solution, rhs );
}

void Richardson::solveFinalize(){
    RichardsonRuntime& runtime = getRuntime();
    if ( runtime.mIterations % 2 ){
        LAMA_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *runtime.mProxyOldSolution = *runtime.mSolution;
    }
    LAMA_LOG_DEBUG( logger, " end solve " )
}

template<typename T>
void Richardson::iterate(){
	typedef T DataType;
    RichardsonRuntime& runtime = getRuntime();

    DataType omega = mOmega.getValue<DataType>();

    const Vector& rhs = *runtime.mRhs;
    const Matrix& A = *runtime.mCoefficients;
    //swap old solution and solution pointer begin
    Vector* ptr_OldSolution = &( *runtime.mProxyOldSolution );
    Vector* ptr_solution = &( *runtime.mSolution );

    runtime.mProxyOldSolution = ptr_solution;
    runtime.mSolution = ptr_OldSolution;

    const Vector& oldSolution = runtime.mProxyOldSolution.getConstReference();

    DenseVector<T> x = A * oldSolution;
   *runtime.mSolution = rhs - x;
	
	if(omega != 1.0)
		*runtime.mSolution = omega * (*runtime.mSolution);

	*runtime.mSolution = oldSolution + (*runtime.mSolution);

    if ( LAMA_LOG_TRACE_ON( logger ) ){
        LAMA_LOG_TRACE( logger, "Solution " << *runtime.mSolution )
        const DenseVector<T>& sol = dynamic_cast<const DenseVector<T>&>( *runtime.mSolution );
        ReadAccess<T> rsol( sol.getLocalValues() );
        std::cout << "Solution: ";
        for ( IndexType i = 0; i < rsol.size(); ++i )
        	std::cout << " " << rsol[i];

        std::cout << std::endl;
    }
}

void Richardson::iterate(){
    switch ( getRuntime().mCoefficients->getValueType() ){
    	case common::scalar::FLOAT:
        	iterate<float>();
        break;
    	case common::scalar::DOUBLE:
       		iterate<double>();
        break;
    	default:
        COMMON_THROWEXCEPTION( "Unsupported ValueType " << getRuntime().mCoefficients->getValueType() )
    }
}

Richardson::RichardsonRuntime& Richardson::getRuntime(){
    return mRichardsonRuntime;
}

const Richardson::RichardsonRuntime& Richardson::getConstRuntime() const{
    return mRichardsonRuntime;
}

SolverPtr Richardson::copy(){
    return SolverPtr( new Richardson( *this ) );
}


} //namespace
