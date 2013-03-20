/**
 * @file SpecialLUSolver.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief SpecialLUSolver.hpp
 * @author Robin Rehrmann
 * @date 18.07.2011
 * $Id$
 */

#ifndef LAMA_LUSOLVER_HPP_
#define LAMA_LUSOLVER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/Solver.hpp>

// others
#include <lama/Communicator.hpp>
#include <lama/DenseMatrix.hpp>

#include <lama/OpenMP/OpenMPBLAS3.hpp>

// tracing
#include <lama/tracing.hpp>

#ifndef LAMA_BUILD_CUDA
typedef int cudaStream_t;
typedef int CUevent;
typedef int CUdevice;
typedef int CUstream;
typedef int CUresult;
typedef int CUDAStreamSyncToken;
typedef std::auto_ptr<CUDAStreamSyncToken> CUDAStreamSyncTokenPtr;
#else
#include <lama/CUDA/CUDAStreamSyncToken.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#endif

namespace lama
{

class LAMA_DLL_IMPORTEXPORT LUSolver: public Solver
{
public:
    /**
     * @brief Creates a solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    LUSolver( const std::string& id );

    /**
     * @brief Create a gaussian solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    LUSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    LUSolver( const LUSolver& other );

    /**
     * @brief LUSolver destructor.
     */
    virtual ~LUSolver();

    /**
     * @brief Used to initialize a gaussian solver with a certain matrix A
     *        from A*u=f.
     *
     * This method initializes a gaussian solver with a certain coefficient-
     * matrix. The given matrix will be overwritten by its lu-factorization.
     *
     * @param coefficients The matrix A from A*u=f.
     */
    virtual void initialize( const Matrix& coefficients );

    void factorMatrixToLU( Matrix& matrix, std::vector<IndexType>& permutation );

    /**
     * @brief Solves the equation system based on the given rhs.
     *
     * The matrix A from A*u=f has to be initialized first!
     * (call gaussianSolver::initialize(matrix) for example). This method uses
     * the lu-factorization, computed in LUSolver::initialize( matrix ),
     * to  calculate  the solution  using forwards  and  backwards application.
     *
     * @param rhs       The right hand side of A*u=f.
     * @param solution  The solution from A*u=f. Mostly used as starting
     *                  solution for an IterativeSolver.
     */
    virtual void solve( Vector& solution, const Vector& rhs );

    void setTileSize( const IndexType tilesize );
    IndexType getTileSize();

    virtual void setDeviceNumber( const IndexType dev );
    virtual IndexType getDeviceNumber();

    struct LUSolverRuntime: SolverRuntime
    {
        LUSolverRuntime();
        virtual ~LUSolverRuntime();

        Matrix* mLUfactorization;
        std::vector<IndexType> mPermutation;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual LUSolverRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const LUSolverRuntime& getConstRuntime() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

protected:

    LUSolverRuntime mLUSolverRuntime;

    IndexType mTilesize;
    IndexType mDev;

    const static double epsilon;

private:
    template<typename T>
    struct lama_swap
    {
        typedef void (*swap_func)( const int, T*, const int, T*, const int );
        ContextPtr ctxt;
        swap_func func;
    };

    template<typename T>
    struct lama_gemm
    {
        typedef void (*gemm_func)(
            const enum CBLAS_ORDER,
            const enum CBLAS_TRANSPOSE,
            const enum CBLAS_TRANSPOSE,
            const int,
            const int,
            const int,
            const T,
            const T*,
            const int,
            const T*,
            const int,
            const T,
            T*,
            const int,
            cudaStream_t );
        typedef void (*__rs)( const CUevent );
        typedef bool (*__qu)( const CUevent );

        cudaStream_t stream;
        gemm_func func;
        __rs record;
        __rs synchronize;
        __qu query;

        static CUDAStreamSyncTokenPtr __syncTok;

        static void __gemm(
            const enum CBLAS_ORDER order,
            const enum CBLAS_TRANSPOSE transa,
            const enum CBLAS_TRANSPOSE transb,
            const int m,
            const int n,
            const int k,
            const T alpha,
            const T* A,
            const int lda,
            const T* B,
            const int ldb,
            const T beta,
            T* C,
            const int ldc,
            cudaStream_t stream );
        static void __recordDef( const CUevent event );
        static bool __queryDef( const CUevent event );
        static void __synchronizeDef( const CUevent event );
        static void __recordCuda( const CUevent event );
        static bool __queryCuda( const CUevent event );
        static void __synchronizeCuda( const CUevent event );
    };

    template<typename T>
    void computeLUFactorization( DenseMatrix<T>& matrix, std::vector<IndexType>& permutation );

    template<typename T>
    void pgetf2(
        const IndexType numBlockRows,
        DenseStorage<T>** const A,
        IndexType* const ipiv,
        const PartitionId ROOT );

    template<typename T>
    void plaswp(
        DenseStorage<T>** const A,
        const PartitionId ROOT,
        const IndexType* const ipiv,
        const IndexType n,
        const lama_swap<T> swap );

    template<typename T>
    IndexType piamax_own(
        const IndexType numBlockCol,
        DenseStorage<T>** const local,
        const IndexType col,
        const IndexType locRow = 0 );

    template<typename T>
    void ptrsm( const enum CBLAS_UPLO uplo, const DenseMatrix<T>& matrix, DenseVector<T>& solution );

    IndexType computeTilesize( IndexType m, IndexType n );

    inline void initializeCommunicator();

    CommunicatorPtr mComm;

    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

// implementation of inner classes.

template<typename T>
CUDAStreamSyncTokenPtr LUSolver::lama_gemm<T>::__syncTok;

template<typename T>
void LUSolver::lama_gemm<T>::__gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_TRANSPOSE transb,
    const int m,
    const int n,
    const int k,
    const T alpha,
    const T* A,
    const int lda,
    const T* B,
    const int ldb,
    const T beta,
    T* C,
    const int ldc,
    cudaStream_t )
{
    LAMA_REGION("GEMM");
    OpenMPBLAS3::gemm( order, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, NULL );
}

template<typename T>
void LUSolver::lama_gemm<T>::__recordDef( const CUevent )
{
}

template<typename T>
void LUSolver::lama_gemm<T>::__synchronizeDef( const CUevent )
{
}

template<typename T>
bool LUSolver::lama_gemm<T>::__queryDef( const CUevent )
{
    return true;
}

template<typename T>
void LUSolver::lama_gemm<T>::__recordCuda( const CUevent event )
{
#ifdef LAMA_BUILD_CUDA
    __syncTok->recordEvent( event );
#endif
}

template<typename T>
void LUSolver::lama_gemm<T>::__synchronizeCuda( const CUevent event )
{
#ifdef LAMA_BUILD_CUDA
    __syncTok->synchronizeEvent( event );
#endif
}

template<typename T>
bool LUSolver::lama_gemm<T>::__queryCuda( const CUevent event )
{
#ifdef LAMA_BUILD_CUDA
    return __syncTok->queryEvent( event );
#endif
    return true;
}

} // namespace LAMA

#endif /* LAMA_LUSOLVER_HPP_ */
