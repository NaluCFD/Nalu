/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ScalarDiffFemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/Hex8FEM.h>

// manage supplemental algorithms that are templated
#include <BuildTemplates.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ScalarDiffFemSuppAlg - hack hex8 FEM heat conduction
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarDiffFemSuppAlg<AlgTraits>::ScalarDiffFemSuppAlg(
  Realm &realm,
  ScalarFieldType *temperature,
  ScalarFieldType *thermalCond)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    temperature_(temperature),
    thermalCond_(thermalCond),
    coordinates_(NULL),
    meFEM_(new Hex8FEM()),
    ipWeight_(&meFEM_->weights_[0])
{
  ThrowRequireMsg(AlgTraits::topo_ == stk::topology::HEX_8,  "FEM_DIFF only available for hexes currently");

  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

  // resize; geometry
  ws_dndx_.resize(AlgTraits::nDim_*AlgTraits::numGp_*AlgTraits::nodesPerElement_);
  ws_deriv_.resize(AlgTraits::nDim_*AlgTraits::numGp_*AlgTraits::nodesPerElement_);
  ws_det_j_.resize(AlgTraits::numGp_);
  ws_shape_function_.resize(AlgTraits::numGp_*AlgTraits::nodesPerElement_);

  // resize; fields
  ws_temperature_.resize(AlgTraits::nodesPerElement_);
  ws_thermalCond_.resize(AlgTraits::nodesPerElement_);
  ws_coordinates_.resize(AlgTraits::nDim_*AlgTraits::nodesPerElement_);
  
  // compute shape function
  meFEM_->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
ScalarDiffFemSuppAlg<AlgTraits>::~ScalarDiffFemSuppAlg()
{
  delete meFEM_;
}

//--------------------------------------------------------------------------
//-------- element_execute -------------------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
void
ScalarDiffFemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  ScratchViews &/*scratchViews*/)
{  
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == AlgTraits::nodesPerElement_ );
  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    
    // gather scalars
    ws_temperature_[ni] = *stk::mesh::field_data(*temperature_, node);
    ws_thermalCond_[ni] = *stk::mesh::field_data(*thermalCond_, node );

    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );

    // gather vectors
    const int offSet = ni*AlgTraits::nDim_;
    for ( int j=0; j < AlgTraits::nDim_; ++j ) {
      ws_coordinates_[offSet+j] = coords[j];
    }
  }

  // compute dndx; det_j comes for free
  double me_error = 0.0;
  meFEM_->grad_op(1, &ws_coordinates_[0], &ws_dndx_[0], &ws_deriv_[0], &ws_det_j_[0], &me_error);

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // compute ip property
    double thermalCondIp = 0.0;
    const int ipNpe = ip*AlgTraits::nodesPerElement_;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = ws_shape_function_[ipNpe+ic];
      thermalCondIp += r*ws_thermalCond_[ic];
    }

    // start the assembly
    const double ipFactor = ws_det_j_[ip]*ipWeight_[ip];

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {

      // offset for row dndx
      const int offSetDnDxIr = AlgTraits::nDim_*AlgTraits::nodesPerElement_*ip + ir*AlgTraits::nDim_;

      // column ic
      double rhsSum = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        // offset for column dndx
        const int offSetDnDxIc = AlgTraits::nDim_*AlgTraits::nodesPerElement_*ip + ic*AlgTraits::nDim_;

        double lhsSum = 0.0;
        double T = ws_temperature_[ic];
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const double fac = ws_dndx_[offSetDnDxIr+j]*thermalCondIp*ws_dndx_[offSetDnDxIc+j];
          lhsSum += fac;
          rhsSum += fac*T;
        }
        lhs[ir*AlgTraits::nodesPerElement_+ic] += lhsSum*ipFactor;
      }
      rhs[ir] -= rhsSum*ipFactor;
    }
  }
}
  
INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ScalarDiffFemSuppAlg);

} // namespace nalu
} // namespace Sierra
