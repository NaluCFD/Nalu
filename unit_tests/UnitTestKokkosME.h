/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTKOKKOSME_H
#define UNITTESTKOKKOSME_H

#include "gtest/gtest.h"
#include "UnitTestUtils.h"

#include "ScratchViews.h"
#include "CopyAndInterleave.h"
#include "ElemDataRequests.h"
#include "AlgTraits.h"
#include "KokkosInterface.h"
#include "SimdInterface.h"

namespace unit_test_utils {

template<typename AlgTraits>
class KokkosMEViews
{
public:
  KokkosMEViews(bool doInit=true, bool doPerturb=false)
    : comm_(MPI_COMM_WORLD),
      meta_(AlgTraits::nDim_),
      bulk_(meta_, comm_)
  {
    if (doInit)
      fill_mesh_and_init_data(doPerturb);
  }

  virtual ~KokkosMEViews() {}

  /** Create a 1-element STK mesh and initialize MasterElement data structures
   */
  void fill_mesh_and_init_data(bool doPerturb=false)
  {
    fill_mesh(doPerturb);
    init_me_data();
  }

  void fill_mesh(bool doPerturb=false)
  {
    if (doPerturb)
      unit_test_utils::create_one_perturbed_element(bulk_, AlgTraits::topo_);
    else
      unit_test_utils::create_one_reference_element(bulk_, AlgTraits::topo_);

    partVec_ = {meta_.get_part("block_1")};
    coordinates_ = static_cast<const VectorFieldType*>(
      meta_.coordinate_field());

    EXPECT_TRUE(coordinates_ != nullptr);
    dataNeeded_.add_coordinates_field(
      *coordinates_, AlgTraits::nDim_, sierra::nalu::CURRENT_COORDINATES);
  }

  void init_me_data()
  {
    // Initialize both surface and volume elements
    meSCS_ = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);
    meSCV_ = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

    // Register them to ElemDataRequests
    dataNeeded_.add_cvfem_surface_me(meSCS_);
    dataNeeded_.add_cvfem_volume_me(meSCV_);

    // Initialize shape function views
    double scs_data[AlgTraits::numScsIp_*AlgTraits::nodesPerElement_];
    meSCS_->shape_fcn(scs_data);
    DoubleType* v_scs_data = &scs_shape_fcn_(0,0);
    for (int i=0; i < (AlgTraits::numScsIp_*AlgTraits::nodesPerElement_); ++i) {
      v_scs_data[i] = scs_data[i];
    }

    double scv_data[AlgTraits::numScvIp_*AlgTraits::nodesPerElement_];
    meSCV_->shape_fcn(scv_data);
    DoubleType* v_scv_data = &scv_shape_fcn_(0,0);
    for (int i=0; i < (AlgTraits::numScvIp_*AlgTraits::nodesPerElement_); ++i) {
      v_scv_data[i] = scv_data[i];
    }
  }

  template<typename LambdaFunction>
  void execute(LambdaFunction func)
  {
    constexpr int simdLen = stk::simd::ndoubles;
    // Do not copy and interleave ME views
    const bool alsoProcessMEViews = false;

    stk::mesh::Selector sel = (
      meta_.locally_owned_part() & stk::mesh::selectUnion(partVec_));

    const int bytes_per_team = 0;
    int bytes_per_thread = sierra::nalu::get_num_bytes_pre_req_data<double>(
      dataNeeded_, AlgTraits::nDim_);
    bytes_per_thread *= 3 * stk::simd::ndoubles;

    const auto& buckets = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);

    auto team_exec = sierra::nalu::get_team_policy(
      buckets.size(), bytes_per_team, bytes_per_thread);

    Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team) {
        auto& b = *buckets[team.league_rank()];
        const auto length = b.size();

        std::vector<std::unique_ptr<sierra::nalu::ScratchViews<double> > > prereqData(simdLen);

        for (int simdIndex=0; simdIndex < simdLen; ++simdIndex) {
          prereqData[simdIndex] = std::unique_ptr<sierra::nalu::ScratchViews<double> >(new sierra::nalu::ScratchViews<double>(team, bulk_, AlgTraits::nodesPerElement_, dataNeeded_));
        }
        sierra::nalu::ScratchViews<DoubleType> simdPrereqData(
          team, bulk_, AlgTraits::nodesPerElement_, dataNeeded_);

        const stk::mesh::Entity* elemNodes[simdLen];
        stk::mesh::Bucket::size_type simdBucketLen = length / simdLen;
        const stk::mesh::Bucket::size_type remainder = length % simdLen;
        if (remainder > 0) simdBucketLen++;

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, simdBucketLen), [&](const size_t& bktIndex){
            int simdElems = simdLen;
            if (length - bktIndex*simdLen < simdLen) {
              simdElems = length - bktIndex*simdLen;
            }

            stk::mesh::Entity element;
            for (int simdIndex=0; simdIndex < simdElems; ++simdIndex) {
              element = b[bktIndex*simdLen + simdIndex];
              elemNodes[simdIndex] = bulk_.begin_nodes(element);
              fill_pre_req_data(dataNeeded_, bulk_, element, *prereqData[simdIndex], alsoProcessMEViews);
            }

            copy_and_interleave(prereqData.data(), simdElems, simdPrereqData,
                                alsoProcessMEViews);
            fill_master_element_views(dataNeeded_, bulk_, simdPrereqData);

            func(simdPrereqData, meSCS_, meSCV_);
          });
      });
  }

  stk::ParallelMachine comm_;
  stk::mesh::MetaData meta_;
  stk::mesh::BulkData bulk_;
  stk::mesh::PartVector partVec_;
  const VectorFieldType* coordinates_{nullptr};

  sierra::nalu::ElemDataRequests dataNeeded_;
  sierra::nalu::MasterElement* meSCV_{nullptr};
  sierra::nalu::MasterElement* meSCS_{nullptr};


  Kokkos::View<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> scv_shape_fcn_ {"scv_shape_function"};
  Kokkos::View<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> scs_shape_fcn_ {"scs_shape_function"};
};

} // namespace unit_test_utils

#endif /* UNITTESTKOKKOSME_H */
