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

const double tol = 1.0e-10;

#ifndef KOKKOS_HAVE_CUDA
//following tests can't run on cuda due to variety of reasons, including
//use of MasterElement functions (defined for host), etc.

void gather_elem_node_coords(int numNodes,
                             int nDim,
                             const stk::mesh::Entity* elemNodes,
                             const VectorFieldType* coordField,
                             double* coordinates)
{
  for(int i=0; i<numNodes; ++i) {
    const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[i]);
    for(int d=0; d<nDim; ++d) {
      coordinates[i*nDim+d] = nodeCoords[d];
    }
  }
}

enum ELEM_DATA_NEEDED {
  NODES = 0,
  COORDS,
  SCS_AREAV,
  SCS_GRAD_OP
};

class ComputedElemData
{
public:
  ComputedElemData(const stk::mesh::BulkData& bulkData,
                   sierra::nalu::MasterElement& meSCS,
                   stk::topology topo,
                   stk::mesh::Entity elem,
                   const std::set<ELEM_DATA_NEEDED>& dataNeeded,
                   const VectorFieldType* coordField)
  {
    int nDim = bulkData.mesh_meta_data().spatial_dimension();
    double scs_error = 0;
    for(ELEM_DATA_NEEDED data : dataNeeded) {
      switch(data)
      {
        case NODES:
           elemNodes = bulkData.begin_nodes(elem);
           haveNodes = true;
           break;
        case COORDS:
           gather_elem_node_coords(topo.num_nodes(), nDim, bulkData.begin_nodes(elem), coordField, coordinates);
           haveCoordinates = true;
           break;
        case SCS_AREAV:
           if (haveCoordinates) {
             meSCS.determinant(1, coordinates, scs_areav, &scs_error);
             haveScsAreav = true;
           }
           break;
        case SCS_GRAD_OP:
           if (haveCoordinates) {
             meSCS.grad_op(1, coordinates, dndx, deriv, det_j, &scs_error);
             haveScsGradOp = true;
           }
           break;
        default: break;
      }
    }
  }

  //the following arrays are sized using 'max' numbers since this is a
  //non-templated version of this class. A templated class could use
  //precisely-sized arrays.
  //
  static const int MAX_NODES_PER_ELEM = 27;
  static const int MAX_NUM_SCS_IP = 216;
  static const int MAX_SPATIAL_DIM = 3;

  bool haveNodes = false;
  const stk::mesh::Entity* elemNodes = nullptr;

  bool haveCoordinates = false;
  double coordinates[MAX_NODES_PER_ELEM * MAX_SPATIAL_DIM];

  bool haveScsAreav = false;
  double scs_areav[MAX_NUM_SCS_IP * MAX_SPATIAL_DIM];

  bool haveScsGradOp = false;
  double dndx[MAX_NODES_PER_ELEM * MAX_NUM_SCS_IP * MAX_SPATIAL_DIM];
  double deriv[MAX_NODES_PER_ELEM * MAX_NUM_SCS_IP * MAX_SPATIAL_DIM];
  double det_j[MAX_NUM_SCS_IP];
};

void element_discrete_laplacian_kernel_3d(
                       sierra::nalu::MasterElement& meSCS,
                       ScalarFieldType* discreteLaplacianOfPressure,
                       ScalarFieldType* nodalPressureField,
                       const ComputedElemData& elemData)
{
    const int nDim = 3;
    const int nodesPerElem = meSCS.nodesPerElement_;
    const int numScsIp = meSCS.numIntPoints_;

    double p_elemNodePressures[nodesPerElem];

    const int* lrscv = meSCS.adjacentNodes();

    const stk::mesh::Entity* elemNodes = elemData.elemNodes;
    for(int n=0; n<nodesPerElem; ++n) {
        const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
        p_elemNodePressures[n] = nodePressure[0];
    }

    const double* scs_areav = elemData.scs_areav;

    const double* p_dndx = elemData.dndx;

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
                    const ComputedElemData& elemData) = 0;
};

class DiscreteLaplacianSuppAlg : public SuppAlg
{
public:
  DiscreteLaplacianSuppAlg(std::set<ELEM_DATA_NEEDED>& dataNeeded,
                           ScalarFieldType* discreteLaplacianOfPressure,
                           ScalarFieldType* nodalPressureField
                          )
   : discreteLaplacianOfPressure_(discreteLaplacianOfPressure),
     nodalPressureField_(nodalPressureField)
  {
    //here's the kinds of element-data we want pre-computed before
    //our elem_execute method is called.
    dataNeeded.insert(NODES);
    dataNeeded.insert(COORDS);
    dataNeeded.insert(SCS_AREAV);
    dataNeeded.insert(SCS_GRAD_OP);
  }

  virtual ~DiscreteLaplacianSuppAlg() {}

  virtual void elem_execute(stk::topology topo,
                    sierra::nalu::MasterElement& meSCS,
                    const ComputedElemData& elemData)
  {
      element_discrete_laplacian_kernel_3d(meSCS,
            discreteLaplacianOfPressure_, nodalPressureField_, elemData);
  }

private:
  ScalarFieldType* discreteLaplacianOfPressure_;
  ScalarFieldType* nodalPressureField_;
};

//=========== Test class that mimics an element alg with supplemental alg ==============
//
class TestElemAlgorithmWithSuppAlg
{
public:
  TestElemAlgorithmWithSuppAlg(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& partVec,
                    const VectorFieldType* coord, ScalarFieldType* discreteLaplacian,
                    ScalarFieldType* nodalPressure)
  : suppAlgs_(), bulkData_(bulk), partVec_(partVec),
    discreteLaplacianOfPressure(discreteLaplacian),
    nodalPressureField(nodalPressure), coordField(coord)
  {}

  void execute()
  {
      const stk::mesh::MetaData& meta = bulkData_.mesh_meta_data();
  
      const stk::mesh::BucketVector& elemBuckets = bulkData_.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  
      kokkos_thread_team_bucket_loop_with_topo(elemBuckets,
          [&](stk::mesh::Entity elem, stk::topology topo, sierra::nalu::MasterElement& meSCS)
          {
             //ComputedElemData must be declared in the inner loop, so that its
             //(compile-time) data is unique for each thread.
             ComputedElemData elemData(bulkData_, meSCS, topo, elem,
                                       dataNeededBySuppAlgs_, coordField);
             for(SuppAlg* alg : suppAlgs_) {
               alg->elem_execute(topo, meSCS, elemData);
             }
          }
      );
  }

  std::vector<SuppAlg*> suppAlgs_;
  std::set<ELEM_DATA_NEEDED> dataNeededBySuppAlgs_;

private:
  stk::mesh::BulkData& bulkData_;
  const stk::mesh::PartVector& partVec_;
  ScalarFieldType* discreteLaplacianOfPressure;
  ScalarFieldType* nodalPressureField;
  const VectorFieldType* coordField;
};


TEST_F(Hex8Mesh, elem_supp_alg)
{
    fill_mesh_and_initialize_test_fields();

    TestElemAlgorithmWithSuppAlg testAlgorithm(bulk, partVec, coordField,
                          discreteLaplacianOfPressure, nodalPressureField);

    //DiscreteLapacianSuppAlg constructor says which data it needs, by inserting
    //things into the 'dataNeededBySuppAlgs_' container.
    SuppAlg* suppAlg = new DiscreteLaplacianSuppAlg(testAlgorithm.dataNeededBySuppAlgs_,
                                           discreteLaplacianOfPressure, nodalPressureField);

    testAlgorithm.suppAlgs_.push_back(suppAlg);

    testAlgorithm.execute();

    check_discrete_laplacian(exactLaplacian);
}

//end of stuff that's ifndef'd for KOKKOS_HAVE_CUDA
#endif

}

