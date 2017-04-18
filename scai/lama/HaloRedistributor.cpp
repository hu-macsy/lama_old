/*
 * HaloRedistributor.cpp
 *
 *  Created on: 18.04.2017
 *      Author: moritzl
 */
#include <vector>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>

#include <scai/dmemo/Redistributor.hpp>

#include "HaloRedistributor.h"

namespace scai {
namespace lama {


HaloRedistributor::HaloRedistributor( DistributionPtr targetDistribution, DistributionPtr sourceDistribution, const Halo &halo) :
		mTargetDistribution(targetDistribution), mSourceDistribution(sourceDistribution), mHalo(halo)
{
	//TODO: maybe put the following sanity check into an assertion, so that it is removed in optimizing builds.
	for (IndexType i = 0; i < targetDistribution->getLocalSize(); i++) {
		IndexType globalI = targetDistribution->local2global(i);
		bool isLocal = sourceDistribution->isLocal(globalI);
		if (!isLocal && halo.global2halo(globalI) == nIndex) {
			throw std::runtime_error("Index " + std::to_string(globalI) + " neither in old distribution, nor in halo.");
		}
	}
}

template<typename T>
void HaloRedistributor::redistributeFromHalo<T>(DenseVector<T>& input, scai::utilskernel::LArray<T>& haloData) {
	SCAI_REGION( "HaloRedistributor.redistributeFromHalo.Vector" )

	scai::dmemo::DistributionPtr mSourceDistribution = input.getDistributionPtr();
	const IndexType newLocalN = mTargetDistribution->getLocalSize();
	LArray<T> newLocalValues;

	{
		scai::hmemo::ReadAccess<T> rOldLocalValues(input.getLocalValues());
		scai::hmemo::ReadAccess<T> rHaloData(haloData);

		scai::hmemo::WriteOnlyAccess<T> wNewLocalValues(newLocalValues, newLocalN);
		for (IndexType i = 0; i < newLocalN; i++) {
			const IndexType globalI = mTargetDistribution->local2global(i);
			if (mSourceDistribution->isLocal(globalI)) {
				wNewLocalValues[i] = rOldLocalValues[mSourceDistribution->global2local(globalI)];
			} else {
				const IndexType localI = mHalo.global2halo(globalI);
				assert(localI != nIndex);
				wNewLocalValues[i] = rHaloData[localI];
			}
		}
	}

	input.swap(newLocalValues, mTargetDistribution);
}

template<typename T>
void HaloRedistributor::redistributeFromHalo(CSRSparseMatrix<T>& matrix, CSRStorage<T>& haloMatrix) {
	SCAI_REGION( "HaloRedistributor.redistributeFromHalo" )

	scai::dmemo::DistributionPtr oldDist = matrix.getRowDistributionPtr();
	//assert that this is the correct distribution

	using scai::utilskernel::LArray;

	const IndexType sourceNumRows = mSourceDistribution->getLocalSize();
	const IndexType targetNumRows = mTargetDistribution->getLocalSize();

	const IndexType globalN = mSourceDistribution->getGlobalSize();
	if (mTargetDistribution->getGlobalSize() != globalN) {
		throw std::runtime_error("Old Distribution has " + std::to_string(globalN) + " values, new distribution has " + std::to_string(mTargetDistribution->getGlobalSize()));
	}

	scai::hmemo::HArray<IndexType> targetIA;
	scai::hmemo::HArray<IndexType> targetJA;
	scai::hmemo::HArray<ValueType> targetValues;

	const CSRStorage<ValueType>& localStorage = matrix.getLocalStorage();

	matrix.setDistributionPtr(mTargetDistribution);

	//check for equality
	if (sourceNumRows == targetNumRows) {
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.equalityCheck" )
		bool allLocal = true;
		for (IndexType i = 0; i < targetNumRows; i++) {
			if (!mSourceDistribution->isLocal(mTargetDistribution->local2global(i))) allLocal = false;
		}
		if (allLocal) {
			//nothing to redistribute and no communication to do either.
			return;
		}
	}

	scai::hmemo::HArray<IndexType> sourceSizes;
	{
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.sourceSizes" )
		scai::hmemo::ReadAccess<IndexType> sourceIA(localStorage.getIA());
		scai::hmemo::WriteOnlyAccess<IndexType> wSourceSizes( sourceSizes, sourceNumRows );
				scai::sparsekernel::OpenMPCSRUtils::offsets2sizes( wSourceSizes.get(), sourceIA.get(), sourceNumRows );
				//allocate
				scai::hmemo::WriteOnlyAccess<IndexType> wTargetIA( targetIA, targetNumRows + 1 );
	}

	scai::hmemo::HArray<IndexType> haloSizes;
	{
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.haloSizes" )
		scai::hmemo::WriteOnlyAccess<IndexType> wHaloSizes( haloSizes, halo.getHaloSize() );
		scai::hmemo::ReadAccess<IndexType> rHaloIA( haloStorage.getIA() );
		scai::sparsekernel::OpenMPCSRUtils::offsets2sizes( wHaloSizes.get(), rHaloIA.get(), halo.getHaloSize() );
	}

	std::vector<IndexType> localTargetIndices;
	std::vector<IndexType> localSourceIndices;
	std::vector<IndexType> localHaloIndices;
	std::vector<IndexType> additionalLocalNodes;
	IndexType numValues = 0;
	{
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.targetIA" )
		scai::hmemo::ReadAccess<IndexType> rSourceSizes(sourceSizes);
		scai::hmemo::ReadAccess<IndexType> rHaloSizes(haloSizes);
				scai::hmemo::WriteAccess<IndexType> wTargetIA( targetIA );

		for (IndexType i = 0; i < targetNumRows; i++) {
			IndexType newGlobalIndex = mTargetDistribution->local2global(i);
			IndexType size;
			if (mSourceDistribution->isLocal(newGlobalIndex)) {
				localTargetIndices.push_back(i);
				const IndexType oldLocalIndex = mSourceDistribution->global2local(newGlobalIndex);
				localSourceIndices.push_back(oldLocalIndex);
				size = rSourceSizes[oldLocalIndex];
			} else {
				additionalLocalNodes.push_back(i);
				const IndexType haloIndex = halo.global2halo(newGlobalIndex);
				assert(haloIndex != nIndex);
				localHaloIndices.push_back(haloIndex);
				size = rHaloSizes[haloIndex];
			}
			wTargetIA[i] = size;
			numValues += size;
		}
		scai::sparsekernel::OpenMPCSRUtils::sizes2offsets( wTargetIA.get(), targetNumRows );

		//allocate
		scai::hmemo::WriteOnlyAccess<IndexType> wTargetJA( targetJA, numValues );
		scai::hmemo::WriteOnlyAccess<ValueType> wTargetValues( targetValues, numValues );
	}

	scai::hmemo::ReadAccess<IndexType> rTargetIA(targetIA);
	assert(rTargetIA.size() == targetNumRows + 1);
	IndexType numLocalIndices = localTargetIndices.size();
	IndexType numHaloIndices = localHaloIndices.size();

	for (IndexType i = 0; i < targetNumRows; i++) {
		assert(rTargetIA[i] <= rTargetIA[i+1]);
				assert(rTargetIA[i] <= numValues);
				//WARNING: the assertion was as below (added '=') but it failed when the last row was empty
				// and rTargetIA[i] = rTargetIA[i+1]
		//assert(rTargetIA[i] < numValues);
	}

	rTargetIA.release();

	{
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.copy" )
		//copying JA array from local matrix and halo
		scai::dmemo::Redistributor::copyV( targetJA, targetIA, LArray<IndexType>(numLocalIndices, localTargetIndices.data()), localStorage.getJA(), localStorage.getIA(), LArray<IndexType>(numLocalIndices, localSourceIndices.data()) );
		scai::dmemo::Redistributor::copyV( targetJA, targetIA, LArray<IndexType>(additionalLocalNodes.size(), additionalLocalNodes.data()), haloStorage.getJA(), haloStorage.getIA(), LArray<IndexType>(numHaloIndices, localHaloIndices.data()) );

		//copying Values array from local matrix and halo
		scai::dmemo::Redistributor::copyV( targetValues, targetIA, LArray<IndexType>(numLocalIndices, localTargetIndices.data()), localStorage.getValues(), localStorage.getIA(), LArray<IndexType>(numLocalIndices, localSourceIndices.data()) );
		scai::dmemo::Redistributor::copyV( targetValues, targetIA, LArray<IndexType>(additionalLocalNodes.size(), additionalLocalNodes.data()), haloStorage.getValues(), haloStorage.getIA(), LArray<IndexType>(numHaloIndices, localHaloIndices.data()) );
	}

	{
		SCAI_REGION( "HaloRedistributor.redistributeFromHalo.setCSRData" )
		//setting CSR data
		matrix.getLocalStorage().setCSRDataSwap(targetNumRows, globalN, numValues, targetIA, targetJA, targetValues, scai::hmemo::ContextPtr());
	}
}

} /* namespace lama */
} /* namespace scai */
