#include <gtest/gtest.h>
#include <limits>
#include <random>

#include "UnitTestUtils.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <Kokkos_Core.hpp>

#include "UnitTestKokkosUtils.h"

namespace {

#ifndef KOKKOS_HAVE_CUDA

void gather_elem_node_field(const stk::mesh::FieldBase& field,
                            stk::mesh::Entity elem,
                            int numNodes,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double*>& shmemView)
{
  unsigned scalarsPerEntity = field.max_size(stk::topology::NODE_RANK);
  for(int i=0; i<numNodes; ++i) {
    const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
    for(unsigned d=0; d<scalarsPerEntity; ++d) {
      shmemView[i*scalarsPerEntity+d] = dataPtr[d];
    }
  }
}

enum ELEM_DATA_NEEDED {
  NODES = 0,
  SCS_AREAV,
  SCS_GRAD_OP
};

struct FieldBaseLess
{
  bool operator()(const stk::mesh::FieldBase* lhs, const stk::mesh::FieldBase* rhs) const
  {
    return lhs->mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal();
  }
};

typedef std::set<const stk::mesh::FieldBase*,FieldBaseLess> FieldSet;

class ElemDataRequests
{
public:
  ElemDataRequests()
  : dataEnums(), fields()
  {
  }

  void add(ELEM_DATA_NEEDED data)
  {
    dataEnums.insert(data);
  }

  void add(const stk::mesh::FieldBase& field)
  {
    fields.insert(&field);
  }

  const std::set<ELEM_DATA_NEEDED>& get_data_enums() const { return dataEnums; }

  const FieldSet& get_fields() const { return fields; }

private:
  std::set<ELEM_DATA_NEEDED> dataEnums;
  FieldSet fields;
};

int get_num_bytes_pre_req_data(int nDim, const ElemDataRequests& dataNeededBySuppAlgs)
{
    const int MAX_NODES_PER_ELEM = 27;
    const int MAX_NUM_SCS_IP = 216;
    int numBytes = 0;

    const FieldSet& neededFields = dataNeededBySuppAlgs.get_fields();
    for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
      stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
      ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
      unsigned scalarsPerEntity = fieldPtr->max_size(fieldEntityRank);
      numBytes += scalarsPerEntity*sizeof(double);
    }

    const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeededBySuppAlgs.get_data_enums();
    for(ELEM_DATA_NEEDED data : dataEnums) {
      switch(data)
      {
        case NODES: numBytes += MAX_NODES_PER_ELEM * sizeof(stk::mesh::Entity);
                    break;
        case SCS_AREAV: numBytes += nDim * MAX_NUM_SCS_IP * sizeof(double);
                    break;
        case SCS_GRAD_OP: numBytes += (MAX_NODES_PER_ELEM*MAX_NUM_SCS_IP*nDim*2 + MAX_NUM_SCS_IP) * sizeof(double);
                    break;
        default: break;
      }
    }

    return numBytes;
}

class ScratchViews
{
public:
  ScratchViews(const TeamHandleType& team,
               const stk::mesh::BulkData& bulkData,
               stk::topology topo,
               const sierra::nalu::MasterElement& meSCS,
               const ElemDataRequests& dataNeeded)
  : elemNodes(), scs_areav(), dndx(), deriv(), det_j()
  {
    int nDim = bulkData.mesh_meta_data().spatial_dimension();
    int nodesPerElem = topo.num_nodes();
    int numScsIp = meSCS.numIntPoints_;

    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    unsigned numRanks = meta.entity_rank_count();
    unsigned numFields = meta.get_fields().size();
    fieldViews.resize(numRanks);
    for(unsigned i=0; i<numRanks; ++i) {
      fieldViews[i].resize(numFields);
    }

    const FieldSet& neededFields = dataNeeded.get_fields();
    for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
      stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
      ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
      unsigned scalarsPerEntity = fieldPtr->max_size(fieldEntityRank);
      unsigned viewLength = fieldEntityRank==stk::topology::ELEM_RANK ? scalarsPerEntity : nodesPerElem*scalarsPerEntity;
      fieldViews[fieldEntityRank][fieldPtr->mesh_meta_data_ordinal()] = get_shmem_view_1D(team, viewLength);
    }

    const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums();
    for(ELEM_DATA_NEEDED data : dataEnums) {
      switch(data)
      {
        case NODES:
           elemNodes = get_entity_shmem_view_1D(team, nodesPerElem);
           break;

        case SCS_AREAV:
           scs_areav = get_shmem_view_1D(team, numScsIp*nDim);
           break;

        case SCS_GRAD_OP:
           dndx = get_shmem_view_1D(team, nodesPerElem*numScsIp*nDim);
           deriv = get_shmem_view_1D(team, nodesPerElem*numScsIp*nDim);
           det_j = get_shmem_view_1D(team, numScsIp);
           break;

        default: break;
      }
    }
  }

  SharedMemView<double*>& get_scratch_view(const stk::mesh::FieldBase& field)
  {
    return fieldViews[field.entity_rank()][field.mesh_meta_data_ordinal()];
  }

  SharedMemView<stk::mesh::Entity*> elemNodes;
  SharedMemView<double*> scs_areav;
  SharedMemView<double*> dndx;
  SharedMemView<double*> deriv;
  SharedMemView<double*> det_j;

  std::vector<std::vector<SharedMemView<double*>>> fieldViews;
};

void fill_pre_req_data(const ElemDataRequests& dataNeeded,
                        const stk::mesh::BulkData& bulkData,
                        stk::topology topo, sierra::nalu::MasterElement& meSCS,
                        stk::mesh::Entity elem,
                        const VectorFieldType* coordField,
                        ScratchViews& prereqData)
{
  double scs_error = 0;
  int nodesPerElem = topo.num_nodes();
  const stk::mesh::Entity* nodes = nullptr;

  const FieldSet& neededFields = dataNeeded.get_fields();
  for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
    SharedMemView<double*>& shmemView = prereqData.get_scratch_view(*fieldPtr);
    if (fieldEntityRank==stk::topology::ELEM_RANK) {
      unsigned len = shmemView.dimension(0);
      double* fieldDataPtr = static_cast<double*>(stk::mesh::field_data(*fieldPtr, elem));
      for(unsigned i=0; i<len; ++i) {
        shmemView(i) = fieldDataPtr[i];
      }
    }
    else if (fieldEntityRank == stk::topology::NODE_RANK) {
      gather_elem_node_field(*fieldPtr, elem, nodesPerElem, bulkData.begin_nodes(elem), shmemView);
    }
    else {
      ThrowRequireMsg(false, "Only node and element fields supported currently.");
    }
  }

  SharedMemView<double*>& coords = prereqData.get_scratch_view(*coordField);
  const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums();
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case NODES:
         nodes = bulkData.begin_nodes(elem);
         for(int i=0; i<nodesPerElem; ++i) {
           prereqData.elemNodes(i) = nodes[i];
         }
         break;

      case SCS_AREAV:
         meSCS.determinant(1, &coords(0), &prereqData.scs_areav(0), &scs_error);
         break;

      case SCS_GRAD_OP:
         meSCS.grad_op(1, &coords(0), &prereqData.dndx(0), &prereqData.deriv(0), &prereqData.det_j(0), &scs_error);
         break;

      default: break;
    }
  }
}

void element_discrete_laplacian_kernel_3d(
                       sierra::nalu::MasterElement& meSCS,
                       const ScalarFieldType* discreteLaplacianOfPressure,
                       const ScalarFieldType* nodalPressureField,
                       ScratchViews& elemData)
{
    const int nDim = 3;
    const int nodesPerElem = meSCS.nodesPerElement_;
    const int numScsIp = meSCS.numIntPoints_;

    const int* lrscv = meSCS.adjacentNodes();

    SharedMemView<double*> nodalPressureView = elemData.get_scratch_view(*nodalPressureField);
    const double* p_elemNodePressures = &nodalPressureView(0);

    const double* scs_areav = &elemData.scs_areav(0);

    const double* p_dndx = &elemData.dndx(0);

    const stk::mesh::Entity* elemNodes = &elemData.elemNodes(0);

    for (int ip = 0; ip < numScsIp; ++ip ) {

      double dpdxIp = 0.0;
      const int ipOffset = nDim*nodesPerElem*ip;
      for ( int ic = 0; ic < nodesPerElem; ++ic) {
        const int offSetDnDx = ipOffset + ic*nDim;
        for ( int j = 0; j < nDim; ++j ) {
          dpdxIp += p_dndx[offSetDnDx+j]*p_elemNodePressures[ic]*scs_areav[ip*nDim+j];
        }
      }
      EXPECT_TRUE(std::abs(dpdxIp) > tol);

      const stk::mesh::Entity lNode = elemNodes[lrscv[2*ip+0]];
      const stk::mesh::Entity rNode = elemNodes[lrscv[2*ip+1]];

      Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, lNode), dpdxIp);
      Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, rNode), -dpdxIp);
    }
}

class SuppAlg
{
public:
  virtual ~SuppAlg(){}

  virtual void elem_execute(stk::topology topo,
                    sierra::nalu::MasterElement& meSCS,
                    ScratchViews& elemData) = 0;
};

class DiscreteLaplacianSuppAlg : public SuppAlg
{
public:
  DiscreteLaplacianSuppAlg(ElemDataRequests& dataNeeded,
                           const VectorFieldType* coordField,
                           const ScalarFieldType* discreteLaplacianOfPressure,
                           const ScalarFieldType* nodalPressureField
                          )
   : discreteLaplacianOfPressure_(discreteLaplacianOfPressure),
     nodalPressureField_(nodalPressureField)
  {
    //here are the element-data pre-requisites we want computed before
    //our elem_execute method is called.
    dataNeeded.add(NODES);
    dataNeeded.add(SCS_AREAV);
    dataNeeded.add(SCS_GRAD_OP);
    dataNeeded.add(*coordField);
    dataNeeded.add(*nodalPressureField);
  }

  virtual ~DiscreteLaplacianSuppAlg() {}

  virtual void elem_execute(stk::topology topo,
                    sierra::nalu::MasterElement& meSCS,
                    ScratchViews& elemData)
  {
      element_discrete_laplacian_kernel_3d(meSCS,
            discreteLaplacianOfPressure_, nodalPressureField_, elemData);
  }

private:
  const ScalarFieldType* discreteLaplacianOfPressure_;
  const ScalarFieldType* nodalPressureField_;
};

//=========== Test class that mimics an element alg with supplemental alg and views ========
//
class TestElemAlgorithmWithSuppAlgViews
{
public:
  TestElemAlgorithmWithSuppAlgViews(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& partVec,
                    const VectorFieldType* coord)
  : suppAlgs_(), bulkData_(bulk), partVec_(partVec),
    coordField(coord)
  {}

  void execute()
  {
      const stk::mesh::MetaData& meta = bulkData_.mesh_meta_data();
  
      const stk::mesh::BucketVector& elemBuckets = bulkData_.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  
      const int bytes_per_team = 0;
      const int bytes_per_thread = get_num_bytes_pre_req_data(meta.spatial_dimension(), dataNeededBySuppAlgs_);
      auto team_exec = get_team_policy(elemBuckets.size(), bytes_per_team, bytes_per_thread);
      Kokkos::parallel_for(team_exec, [&](const TeamHandleType& team)
      {
          const stk::mesh::Bucket& bkt = *elemBuckets[team.league_rank()];
          stk::topology topo = bkt.topology();
          sierra::nalu::MasterElement& meSCS = *unit_test_utils::get_surface_master_element(topo);

          ScratchViews prereqData(team, bulkData_, topo, meSCS, dataNeededBySuppAlgs_);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& jj)
          {
             fill_pre_req_data(dataNeededBySuppAlgs_, bulkData_, topo, meSCS,
                                bkt[jj], coordField, prereqData);

             for(SuppAlg* alg : suppAlgs_) {
               alg->elem_execute(topo, meSCS, prereqData);
             }
          });
      });
  }

  std::vector<SuppAlg*> suppAlgs_;
  ElemDataRequests dataNeededBySuppAlgs_;

private:
  stk::mesh::BulkData& bulkData_;
  const stk::mesh::PartVector& partVec_;
  const VectorFieldType* coordField;
};


TEST_F(Hex8Mesh, elem_supp_alg_views)
{
    fill_mesh_and_initialize_test_fields();

    TestElemAlgorithmWithSuppAlgViews testAlgorithm(bulk, partVec, coordField);

    //DiscreteLapacianSuppAlg constructor says which data it needs, by inserting
    //things into the 'dataNeededBySuppAlgs_' container.
    SuppAlg* suppAlg = new DiscreteLaplacianSuppAlg(testAlgorithm.dataNeededBySuppAlgs_,
                                           coordField,
                                           discreteLaplacianOfPressure, nodalPressureField);

    testAlgorithm.suppAlgs_.push_back(suppAlg);

    testAlgorithm.execute();

    check_discrete_laplacian(exactLaplacian);

    delete suppAlg;
}

//end of stuff that's ifndef'd for KOKKOS_HAVE_CUDA
#endif

}

