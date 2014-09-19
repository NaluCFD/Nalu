/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleHeatCondWallSolverAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleHeatCondWallSolverAlgorithm - provides h(T-Too) flux
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleHeatCondWallSolverAlgorithm::AssembleHeatCondWallSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  bool useShifted)
  : SolverAlgorithm(realm, part, eqSystem),
    useShifted_(useShifted)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  heatTransCoeff_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_transfer_coefficient");
  referenceTemp_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "reference_temperature");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleHeatCondWallSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildFaceToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleHeatCondWallSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<stk::mesh::Entity> connected_nodes;

  // nodal fields to gather
  std::vector<double> ws_htc;
  std::vector<double> ws_reference_temperature;
  std::vector<double> ws_temperature;

  // geometry related to populate
  std::vector<double> ws_shape_function;

  // setup for buckets; union parts and ask for locally owned
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract master element specifics
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsIp = meFC->numIntPoints_;

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nodesPerFace;
    const int rhsSize = nodesPerFace;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    connected_nodes.resize(nodesPerFace);

    // algorithm related
    ws_htc.resize(nodesPerFace);
    ws_reference_temperature.resize(nodesPerFace);
    ws_temperature.resize(nodesPerFace);
    ws_shape_function.resize(numScsIp*nodesPerFace);

    // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_htc = &ws_htc[0];
    double *p_reference_temperature = &ws_reference_temperature[0];
    double *p_temperature = &ws_temperature[0];
    double *p_shape_function = &ws_shape_function[0];

    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_shape_function[0]);
    else
      meFC->shape_fcn(&p_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // face data
      double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      // face node relations for nodal gather
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);
      
      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {

        // get the node and form connected_node
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;

        // gather scalar
        p_htc[ni] = *stk::mesh::field_data(*heatTransCoeff_, node);
        p_reference_temperature[ni] = *stk::mesh::field_data(*referenceTemp_, node);
        p_temperature[ni] = *stk::mesh::field_data(*temperature_, node);
      }

      // start the assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {
	
        double magA = 0.0;
        for ( int j=0; j < nDim; ++j ) {
          magA += areaVec[ip*nDim+j]*areaVec[ip*nDim+j];
        }
        magA = std::sqrt(magA);
	
        const int nn = ip;
        const int offSet = ip*nodesPerFace;
	
        // form boundary ip values
        double htcBip = 0.0;
        double tRefBip = 0.0;
        double tBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_shape_function[offSet+ic];
          htcBip += r*p_htc[ic];
          tRefBip += r*p_reference_temperature[ic];
          tBip += r*p_temperature[ic];
        }

        // form convection and rhs contribution
        const double convection = htcBip*(tRefBip - tBip)*magA;
        p_rhs[nn] += convection;
	
        // sensitivities
        const int rowR = nn*nodesPerFace;
        const double lhsFac = htcBip*magA;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_shape_function[offSet+ic];
          p_lhs[rowR+ic] += r*lhsFac;
        }
       
      }
      
      apply_coeff(connected_nodes, rhs, lhs, __FILE__);

    }
  }
}

} // namespace nalu
} // namespace Sierra
