#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_simd/Simd.hpp>
#include <KokkosInterface.h>

#include "UnitTestUtils.h"
#include <stk_mesh/base/GetEntities.hpp>

#include <SolutionOptions.h>
#include <nso/MomentumNSOElemKernel.h>
#include <ElemDataRequests.h>
#include <ScratchViews.h>
#include <master_element/MasterElement.h>

#include <limits>
#include <vector>

template<typename T>
using AlignedVector = std::vector<T, non_std::AlignedAllocator<T,64> >;



void initialize_coords(Kokkos::View<double**>& coords)
{
    coords(0,0) = -0.5; coords(0,1) = -0.5;  coords(0,2) = 0.5;
    coords(1,0) = -0.5; coords(1,1) = -0.5; coords(1,2) = -0.5;
    coords(2,0) = -0.5; coords(2,1) = 0.5; coords(2,2) = -0.5;
    coords(3,0) = -0.5; coords(3,1) = 0.5; coords(3,2) = 0.5;
    coords(4,0) = 0.5; coords(4,1) = -0.5; coords(4,2) = 0.5;
    coords(5,0) = 0.5; coords(5,1) = -0.5; coords(5,2) = -0.5;
    coords(6,0) = 0.5; coords(6,1) = 0.5; coords(6,2) = -0.5;
    coords(7,0) = 0.5; coords(7,1) = 0.5; coords(7,2) = 0.5;
}

void zero(Kokkos::View<double**>& lhs, Kokkos::View<double*>& rhs)
{
  for(size_t i=0; i<lhs.dimension(0); ++i) {
    for(size_t j=0; j<lhs.dimension(1); ++j) {
      lhs(i,j) = 0.0;
    }
  }
  for(size_t i=0; i<rhs.dimension(0); ++i) {
    rhs(i) = 0.0;
  }
}

void zero(sierra::nalu::SharedMemView<DoubleType**>& lhs, sierra::nalu::SharedMemView<DoubleType*>& rhs)
{
  for(size_t i=0; i<lhs.dimension(0); ++i) {
    for(size_t j=0; j<lhs.dimension(1); ++j) {
      lhs(i,j) = 0.0;
    }
  }
  for(size_t i=0; i<rhs.dimension(0); ++i) {
    rhs(i) = 0.0;
  }
}

void zero(Kokkos::View<stk::simd::Double**>& lhs, Kokkos::View<stk::simd::Double*>& rhs)
{
  for(size_t i=0; i<lhs.dimension(0); ++i) {
    for(size_t j=0; j<lhs.dimension(1); ++j) {
      lhs(i,j) = 0.0;
    }
  }
  for(size_t i=0; i<rhs.dimension(0); ++i) {
    rhs(i) = 0.0;
  }
}


//TEST_F(Hex8MeshWithNSOFields, twoMomentumKernels)
//{
//  fill_mesh_and_initialize_test_fields("generated:50x50x50");
//
//  const VectorFieldType* velocityNp1 = &velocity->field_of_state(stk::mesh::StateNP1);
//  const ScalarFieldType* densityNp1 = &density->field_of_state(stk::mesh::StateNP1);
//  sierra::nalu::HexSCS hex8SCS;
//  int integPts = hex8SCS.numIntPoints_;
//  int nodesPerElem = hex8SCS.nodesPerElement_;
//  int spatialDim = hex8SCS.nDim_;
//
//  stk::mesh::Selector all = meta.universal_part();
//  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, all);
//  const int numElems = stk::mesh::count_selected_entities(all, elemBuckets);
//
//  int simdLen = 1;
//  Data data; 
//  data.setup(integPts, nodesPerElem, spatialDim, spatialDim, simdLen, &hex8SCS);
//  DataSimd<double> dataSimd(data);
//
//  double startTime = stk::wall_time();
//  unsigned calls = 0;
//  for(const stk::mesh::Bucket* bptr : elemBuckets) {
//    const stk::mesh::Bucket& bkt = *bptr;
//    for(size_t i=0; i<bkt.size(); ++i) {
//      stk::mesh::Entity elem = bkt[i];
//      data.gather_and_compute(0, elem, coordField, scalarQ, diffFluxCoeff, massFlowRate, viscosity,
//                           velocityNp1, densityNp1, pressure);
//      dataSimd.copy_and_interleave(data);
//      zero(dataSimd.lhs, dataSimd.rhs);
//
//      ++calls;
//      momentumAdvDiffElem(integPts, nodesPerElem, spatialDim, dataSimd);
//      momentumNSOElem(integPts, nodesPerElem, spatialDim, dataSimd);
//    }
//  }
//
//  double elapsedTimeNoSimd = stk::wall_time() - startTime;
//  std::cout<<"numElems: "<<numElems<<", elapsedTime Hex8MeshWithNSOFields.twoMomentumKernels: "<<elapsedTimeNoSimd<<", calls: "<<calls<<std::endl;
//}

TEST_F(Hex8MeshWithNSOFields, twoMomentumKernelsSimd)
{
  fill_mesh_and_initialize_test_fields("generated:20x20x20");

  sierra::nalu::SolutionOptions solnOpts;
  solnOpts.meshMotion_ = false;
  solnOpts.meshDeformation_ = false;
  solnOpts.externalMeshDeformation_ = false;
  solnOpts.includeDivU_ = 0.0;

  int spatialDim = bulk.mesh_meta_data().spatial_dimension();

  sierra::nalu::HexSCS hex8SCS;

  sierra::nalu::ElemDataRequests dataNeeded;
  dataNeeded.add_coordinates_field(*coordField, spatialDim, sierra::nalu::CURRENT_COORDINATES);
  dataNeeded.add_cvfem_surface_me(&hex8SCS);

  sierra::nalu::MomentumNSOElemKernel<sierra::nalu::AlgTraitsHex8> nsoAlg(
     bulk, solnOpts, velocity, Gju, viscosity, 0.0, 0.0, dataNeeded);

  int nodesPerElem = hex8SCS.nodesPerElement_;

  const int rhsSize = nodesPerElem*spatialDim;
  const int lhsSize = rhsSize * rhsSize;
  const int scratchIdsSize = rhsSize;

  const int simdLen = stk::simd::ndoubles; //only turn this on if you also switch DoubleType

  const int bytes_per_team = 0;
  int bytes_per_thread = (rhsSize + lhsSize) * sizeof(double) + scratchIdsSize*sizeof(int) +
     sierra::nalu::get_num_bytes_pre_req_data(dataNeeded, spatialDim);
  bytes_per_thread *= 2*simdLen;
  

  stk::mesh::Selector all = meta.universal_part();
  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, all);
  const int numElems = stk::mesh::count_selected_entities(all, elemBuckets);

  auto team_exec = sierra::nalu::get_team_policy(elemBuckets.size(), bytes_per_team, bytes_per_thread);

  std::vector<sierra::nalu::ScratchViews<double>*> prereqData(simdLen, nullptr);

  double startTime = stk::wall_time();
  unsigned calls = 0;
  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team)
  {
    const stk::mesh::Bucket &bkt = *elemBuckets[team.league_rank()];

    for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
      delete prereqData[simdIndex];
      prereqData[simdIndex] = new sierra::nalu::ScratchViews<double>(team, bulk, bkt.topology(), dataNeeded);
    }

    sierra::nalu::ScratchViews<DoubleType> prereqDataSimd(team, bulk, bkt.topology(), dataNeeded);
    sierra::nalu::SharedMemView<int*> scratchIds  = sierra::nalu::get_int_shmem_view_1D(team, scratchIdsSize);
    sierra::nalu::SharedMemView<DoubleType*> rhs  = sierra::nalu::get_shmem_view_1D<DoubleType>(team, rhsSize);
    sierra::nalu::SharedMemView<DoubleType**> lhs = sierra::nalu::get_shmem_view_2D<DoubleType>(team, rhsSize, rhsSize);

    for(size_t bktIndex=0; bktIndex<bkt.size(); bktIndex += simdLen)
    {
      int nelem = simdLen;
      if (bkt.size() - bktIndex < simdLen) {
        nelem = bkt.size() - bktIndex;
      }
      if (nelem < 0 || nelem > simdLen) {
        std::cout<<" nelem == "<<nelem<<" shouldn't happen!!!"<<std::endl;
        break;
      }
 
      for(int simdElemIndex=0; simdElemIndex<nelem; ++simdElemIndex) {
        stk::mesh::Entity elem = bkt[bktIndex+simdElemIndex];
        fill_pre_req_data(dataNeeded, bulk, bkt.topology(), elem, *prereqData[simdElemIndex]);
      }
      copy_and_interleave(prereqData, nelem, prereqDataSimd);
      zero(lhs, rhs);

      ++calls;
      nsoAlg.execute(lhs, rhs, prereqDataSimd);

//      momentumAdvDiffElem(integPts, nodesPerElem, spatialDim, prereqDataSimd);
//      momentumNSOElem(integPts, nodesPerElem, spatialDim, prereqDataSimd);
    }
  });

  double elapsedTimeSimd = stk::wall_time() - startTime;
  std::cout<<"numElems: "<<numElems<<", elapsedTime Hex8MeshWithNSOFields.twoMomentumKernels: "<<elapsedTimeSimd<<", calls: "<<calls<<std::endl;
}

