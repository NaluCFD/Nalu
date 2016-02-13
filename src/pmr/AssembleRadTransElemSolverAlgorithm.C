/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <pmr/AssembleRadTransElemSolverAlgorithm.h>
#include <pmr/RadiativeTransportEquationSystem.h>
#include <SolverAlgorithm.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleRadTransElemSolverAlgorithm - add LHS/RHS for Intensity (SUCV)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleRadTransElemSolverAlgorithm::AssembleRadTransElemSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  RadiativeTransportEquationSystem *radEqSystem)
  : SolverAlgorithm(realm, part, radEqSystem),
    radEqSystem_(radEqSystem),
    intensity_(NULL),
    coordinates_(NULL),
    absorption_(NULL),
    scattering_(NULL),
    scalarFlux_(NULL),
    radiationSource_(NULL),
    dualNodalVolume_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  absorption_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "absorption_coefficient");
  scattering_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "scattering_coefficient");
  scalarFlux_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "scalar_flux");
  radiationSource_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "radiation_source");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
 
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleRadTransElemSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleRadTransElemSolverAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // use edge-based length scale
  const bool useEdgeH = true;

  // extract current ordinate direction
  std::vector<double> Sk(nDim,0.0);
  radEqSystem_->get_current_ordinate(&Sk[0]);
  const double *p_Sk = &Sk[0];
  intensity_ = radEqSystem_->get_intensity();
  
  const double invPi = 1.0/(std::acos(-1.0));

   // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_intensity;
  std::vector<double> ws_absorption;
  std::vector<double> ws_scattering;
  std::vector<double> ws_scalarFlux;
  std::vector<double> ws_radiationSource;
  std::vector<double> ws_dualVolume;

  // geometry related to populate
  std::vector<double> ws_scs_areav;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;
  std::vector<double> ws_shape_function;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const size_t length   = b.size();

    // extract master element
    MasterElement *meSCS = realm_.get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();

    // resize some things; matrix related
    const int lhsSize = nodesPerElement*nodesPerElement;
    const int rhsSize = nodesPerElement;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerElement);

    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_intensity.resize(nodesPerElement);
    ws_absorption.resize(nodesPerElement);
    ws_scattering.resize(nodesPerElement);
    ws_scalarFlux.resize(nodesPerElement);
    ws_radiationSource.resize(nodesPerElement);
    ws_dualVolume.resize(nodesPerElement);
    ws_scs_areav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_det_j.resize(numScsIp);
    ws_shape_function.resize(numScsIp*nodesPerElement);

     // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_intensity = &ws_intensity[0];
    double *p_absorption = &ws_absorption[0];
    double *p_scattering = &ws_scattering[0];
    double *p_scalarFlux = &ws_scalarFlux[0];
    double *p_radiationSource = &ws_radiationSource[0];
    double *p_dualVolume = &ws_dualVolume[0];
    double *p_scs_areav = &ws_scs_areav[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];

    meSCS->shape_fcn(&p_shape_function[0]);

    for ( size_t k = 0 ; k < length ; ++k ) {

        // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // get elem and its node relations
      unsigned elem_offset = k;

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(elem_offset);
      int num_nodes = b.num_nodes(elem_offset);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // set connected nodes
        connected_nodes[ni] = node;

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates_, node);

        // gather scalars
        p_intensity[ni]   = *stk::mesh::field_data(*intensity_, node);
        p_absorption[ni]  = *stk::mesh::field_data(*absorption_, node );
        p_scattering[ni]  = *stk::mesh::field_data(*scattering_, node );
        p_scalarFlux[ni]  = *stk::mesh::field_data(*scalarFlux_, node );
        p_radiationSource[ni] = *stk::mesh::field_data(*radiationSource_, node );
        p_dualVolume[ni]  = *stk::mesh::field_data(*dualNodalVolume_, node );

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // compute geometry
      double scs_error = 0.0;
      meSCS->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

      // compute dndx
      meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);

      for ( int ip = 0; ip < numScsIp; ++ip ) {
	
        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];

        // corresponding matrix rows
        int rowL = il*nodesPerElement;
        int rowR = ir*nodesPerElement;

        // form sj*njdS (part of the lhs for central term; I*sj*njdS)
        double sjaj = 0.0;
        double asq = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double aj = p_scs_areav[ip*nDim+j];
          sjaj += p_Sk[j]*aj;
          asq += aj*aj;
        }
        const double aMag = std::sqrt(asq);

        // integration point interpolation
        double Iscs = 0.0;
        double extCoeffscs = 0.0;
        double ePscs = 0.0;
        double isotropicScatterscs = 0.0;
        double dualNodalVscs = 0.0;
        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSet+ic];
          
          // save of some variables
          const double I = p_intensity[ic];
          const double mua = p_absorption[ic];
          const double mus = p_scattering[ic];

          // interpolation to scs
          Iscs += r*I;
          extCoeffscs += r*(mua+mus);
          ePscs += r*p_radiationSource[ic];
          isotropicScatterscs += r*mus*p_scalarFlux[ic]/4.0*invPi;
          dualNodalVscs += r*p_dualVolume[ic];

          // assemble I*sj*njdS to lhs; left/right
          p_lhs[rowL+ic] += sjaj*r;
          p_lhs[rowR+ic] -= sjaj*r;
        }

        // rhs residual for I*sj*njdS
        p_rhs[il] -= Iscs*sjaj;
        p_rhs[ir] += Iscs*sjaj;

        // now work on SUCV stabilization terms; needed tau, hence second ic loop
        double h_edge = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double nj = p_scs_areav[ip*nDim+j]/aMag;
          const double dxj = p_coordinates[ir*nDim+j]-p_coordinates[il*nDim+j];
          h_edge += nj*dxj;
        }

        // alternative h
        const double h_vol = std::pow(dualNodalVscs, 1.0/(double)nDim);

        // form tau
        const double h = (useEdgeH) ? h_edge : h_vol;
        const double tau = std::sqrt(1.0/((2.0/h)*(2.0/h) + extCoeffscs*extCoeffscs));
	
        double sidIdxi = 0.0;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const double r = p_shape_function[offSet+ic];

          // save of some variables
          const double I = p_intensity[ic];
          
          // SUCV -tau*sj*aj*(mua+mus)*I term; left/right (residual below)
          p_lhs[rowL+ic] += -tau*sjaj*r*extCoeffscs;
          p_lhs[rowR+ic] -= -tau*sjaj*r*extCoeffscs;
	  
          // SUCV diffusion-like term; -tau*si*dI/dxi*sjaj (residual below)
          double lhsfac = 0.0;
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            const double sjdNj = p_Sk[j]*p_dndx[offSetDnDx+j];
            sidIdxi += sjdNj*I;
            lhsfac += -sjdNj;
          }
          p_lhs[rowL+ic] += tau*sjaj*lhsfac;
          p_lhs[rowR+ic] -= tau*sjaj*lhsfac;
	  
        }
	
        // full sucv residual
	const double residual = -tau*sjaj*(sidIdxi + extCoeffscs*Iscs - ePscs - isotropicScatterscs);
	
        // residual; left and right
        p_rhs[il] -= residual;
        p_rhs[ir] += residual;
	
      }

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
    }

  }
}

} // namespace nalu
} // namespace Sierra
