/*
 * HaloRedistributor.h
 *
 *  Created on: 18.04.2017
 *      Author: moritzl
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Halo.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/hmemo/HArray.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/unique_ptr.hpp>

namespace scai {

using dmemo::DistributionPtr;
using dmemo::Halo;
using utilskernel::LArray;

namespace lama {

class COMMON_DLL_IMPORTEXPORT HaloRedistributor: public scai::common::Printable
{
public:
	HaloRedistributor( DistributionPtr targetDistribution, DistributionPtr sourceDistribution, const Halo &halo);
	virtual ~HaloRedistributor();

    template<typename T>
	void redistributeFromHalo(CSRSparseMatrix<T>& matrix, CSRStorage<T>& haloMatrix);

	template<typename T>
	void redistributeFromHalo(DenseVector<T>& input, scai::utilskernel::LArray<T>& haloData);


private:
	DistributionPtr mTargetDistribution;
	DistributionPtr mSourceDistribution;

	Halo mHalo;

};

} /* namespace lama */
} /* namespace scai */
