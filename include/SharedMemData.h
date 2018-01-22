/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SharedMemData_h
#define SharedMemData_h

#include <stk_mesh/base/BulkData.hpp>

#include <ElemDataRequests.h>
#include <KokkosInterface.h>

#include <memory>

namespace sierra{
namespace nalu{

static constexpr int simdLen = stk::simd::ndoubles;

struct SharedMemData {
    SharedMemData(const sierra::nalu::TeamHandleType& team,
         const stk::mesh::BulkData& bulk,
         const ElemDataRequests& dataNeededByKernels,
         unsigned nodesPerEntity,
         unsigned rhsSize)
     : simdPrereqData(team, bulk, nodesPerEntity, dataNeededByKernels)
    {
        for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
          prereqData[simdIndex] = std::unique_ptr<ScratchViews<double> >(new ScratchViews<double>(team, bulk, nodesPerEntity, dataNeededByKernels));
        }
        simdrhs = get_shmem_view_1D<DoubleType>(team, rhsSize);
        simdlhs = get_shmem_view_2D<DoubleType>(team, rhsSize, rhsSize);
        rhs = get_shmem_view_1D<double>(team, rhsSize);
        lhs = get_shmem_view_2D<double>(team, rhsSize, rhsSize);

        scratchIds = get_int_shmem_view_1D(team, rhsSize);
        sortPermutation = get_int_shmem_view_1D(team, rhsSize);
    }

    int numSimdElems;
    std::unique_ptr<ScratchViews<double>> prereqData[simdLen];
    ScratchViews<DoubleType> simdPrereqData;
    SharedMemView<DoubleType*> simdrhs;
    SharedMemView<DoubleType**> simdlhs;
    SharedMemView<double*> rhs;
    SharedMemView<double**> lhs;

    SharedMemView<int*> scratchIds;
    SharedMemView<int*> sortPermutation;
};

} // namespace nalu
} // namespace Sierra

#endif
