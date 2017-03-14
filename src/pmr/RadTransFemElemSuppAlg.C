/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <pmr/RadTransFemElemSuppAlg.h>
#include <pmr/RadiativeTransportEquationSystem.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/Hex8FEM.h>

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
// RadTransFemElemSuppAlg - hack hex8 FEM SUPG
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
RadTransFemElemSuppAlg::RadTransFemElemSuppAlg(
  Realm &realm,
  RadiativeTransportEquationSystem *radEqSystem)
  : SupplementalAlgorithm(realm),
    radEqSystem_(radEqSystem),
    bulkData_(&realm.bulk_data()),
    intensity_(NULL),
    absorption_(NULL),
    scattering_(NULL),
    scalarFlux_(NULL),
    radiationSource_(NULL),
    coordinates_(NULL),
    meFEM_(new Hex8FEM()),
    ipWeight_(&meFEM_->weights_[0]),
    invPi_(1.0/std::acos(-1.0)),
    nDim_(realm.meta_data().spatial_dimension()),
    rowSumLump_(0.0),
    consistentMass_(1.0-rowSumLump_),
    linearNorm_(true),
    useUpper_(true)
{
  // save off fields; for non-BDF2 gather in state N for Nm1 (gamma3_ will be zero)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  absorption_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "absorption_coefficient");
  scattering_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "scattering_coefficient");
  scalarFlux_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "scalar_flux");
  radiationSource_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "radiation_source");
 
  // fixed size
  ws_Sk_.resize(nDim_);
  
  // tell the user what norm
  if ( linearNorm_ ) 
    NaluEnv::self().naluOutputP0() << "Linear Norm Activated " ;
  else 
    NaluEnv::self().naluOutputP0() << "Euclidean Norm Activated " ;

  // and g^ij or g_ij
  if ( useUpper_ )
    NaluEnv::self().naluOutputP0() << "with g^ij form " << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "with g_ij form " << std::endl;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
RadTransFemElemSuppAlg::~RadTransFemElemSuppAlg()
{
  delete meFEM_;
}

//--------------------------------------------------------------------------
//-------- elem_resize -----------------------------------------------------
//--------------------------------------------------------------------------
void
RadTransFemElemSuppAlg::elem_resize(
  MasterElement */*meFEM_*/,
  MasterElement */*meSCV*/)
{
  const int nodesPerElement = meFEM_->nodesPerElement_;
  const int numIp = meFEM_->numIntPoints_;

  // resize; geometry
  ws_dndx_.resize(nDim_*numIp*nodesPerElement);
  ws_deriv_.resize(nDim_*numIp*nodesPerElement);
  ws_det_j_.resize(numIp);
  ws_shape_function_.resize(numIp*nodesPerElement);
  ws_gUpper_.resize(nDim_*nDim_*numIp); // g^ij (covariant)
  ws_gLower_.resize(nDim_*nDim_*numIp); // g_ij (contravariat)

  // resize; fields
  ws_intensity_.resize(nodesPerElement);
  ws_absorption_.resize(nodesPerElement);
  ws_scattering_.resize(nodesPerElement);
  ws_scalarFlux_.resize(nodesPerElement);
  ws_radiationSource_.resize(nodesPerElement);
  ws_coordinates_.resize(nDim_*nodesPerElement);
  
  // compute shape function
  meFEM_->shape_fcn(&ws_shape_function_[0]);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
RadTransFemElemSuppAlg::setup()
{
  // extract current ordinate direction
  radEqSystem_->get_current_ordinate(&ws_Sk_[0]);
  intensity_ = radEqSystem_->get_intensity();
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
RadTransFemElemSuppAlg::elem_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity element,
  MasterElement */*meSCS_*/,
  MasterElement */*meSCV*/)
{
  // details on this element topo
  const int nodesPerElement = meFEM_->nodesPerElement_;
  const int numIp = meFEM_->numIntPoints_;
  
  // gather
  stk::mesh::Entity const *  node_rels = bulkData_->begin_nodes(element);
  int num_nodes = bulkData_->num_nodes(element);

  // sanity check on num nodes
  ThrowAssert( num_nodes == nodesPerElement );

  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = node_rels[ni];
    
    // gather scalars
    ws_intensity_[ni]   = *stk::mesh::field_data(*intensity_, node);
    ws_absorption_[ni]  = *stk::mesh::field_data(*absorption_, node );
    ws_scattering_[ni]  = *stk::mesh::field_data(*scattering_, node );
    ws_scalarFlux_[ni]  = *stk::mesh::field_data(*scalarFlux_, node );
    ws_radiationSource_[ni] = *stk::mesh::field_data(*radiationSource_, node );

    // pointers to real data
    const double * coords = stk::mesh::field_data(*coordinates_, node );

    // gather vectors
    const int offSet = ni*nDim_;
    for ( int j=0; j < nDim_; ++j ) {
      ws_coordinates_[offSet+j] = coords[j];
    }
  }

  // compute dndx and gij
  double me_error = 0.0; 
  meFEM_->grad_op(1, &ws_coordinates_[0], &ws_dndx_[0], &ws_deriv_[0], &ws_det_j_[0], &me_error);
  meFEM_->gij(&ws_coordinates_[0], &ws_gUpper_[0], &ws_gLower_[0], &ws_deriv_[0]);

  for ( int ip = 0; ip < numIp; ++ip ) {

    // pointer to glowerij
    const double *p_gUpper = &ws_gUpper_[nDim_*nDim_*ip];
    const double *p_gLower = &ws_gLower_[nDim_*nDim_*ip];
   
    double Iip = 0.0;
    double extCoeffIp = 0.0;
    double radSrcIp = 0.0;
    double isotropicScatterIp = 0.0;
    double sdIdx = 0.0;
    const int ipNpe = ip*nodesPerElement;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      const double r = ws_shape_function_[ipNpe+ic];

      // save of some variables
      const double I = ws_intensity_[ic];
      const double mua = ws_absorption_[ic];
      const double mus = ws_scattering_[ic];

      // interpolation to ip
      Iip += r*I;
      extCoeffIp += r*(mua+mus);
      radSrcIp += r*ws_radiationSource_[ic];
      isotropicScatterIp += r*mus*ws_scalarFlux_[ic]/4.0*invPi_;

      // compute ip derivatives
      const int offSetDnDx = nDim_*nodesPerElement*ip + ic*nDim_;
      for ( int j = 0; j < nDim_; ++j ) {
        const double dnj = ws_dndx_[offSetDnDx+j];
        sdIdx += I*dnj*ws_Sk_[j];
      }
    }

    // compute "flow-aligned" length scale
    double flowAlignedUp = 0.0;
    double flowAlignedLow = 0.0;
    for ( int i = 0; i < nDim_; ++i ) {
      for ( int j = 0; j < nDim_; ++j ) {
        flowAlignedUp  += ws_Sk_[i]*p_gUpper[i*nDim_+j]*ws_Sk_[j];
        flowAlignedLow  += ws_Sk_[i]*p_gLower[i*nDim_+j]*ws_Sk_[j];
      }
    }      

    // determine inverse length scale; ad-hoc works, however, not formally correct
    double tau = 0.0;
    if ( useUpper_ ) {
      // g^ij form; contravariant
      if ( linearNorm_ )
        tau = 1.0/(1.0/std::sqrt(flowAlignedUp) + extCoeffIp);
      else
        tau = 1.0/(std::sqrt(1.0/flowAlignedUp + extCoeffIp*extCoeffIp));
    }
    else {
      // g_ij form; covariant
      if ( linearNorm_ )
        tau = 1.0/(std::sqrt(flowAlignedLow) + extCoeffIp);
      else
        tau = 1.0/(std::sqrt(flowAlignedLow + extCoeffIp*extCoeffIp));
    }
      
    // compute scaled residual for this ip
    const double residual = sdIdx + extCoeffIp*Iip - radSrcIp - isotropicScatterIp;

    // start the assembly
    const double ipFactor = ws_det_j_[ip] * ipWeight_[ip];

    // row ir
    for ( int ir = 0; ir < nodesPerElement; ++ir) {

      // offset for row dndx
      const int offSetDnDxIr = nDim_*nodesPerElement*ip + ir*nDim_;
      const double rIr = ws_shape_function_[ipNpe+ir];

      double sdWdx = 0.0;
      for ( int j = 0; j < nDim_; ++j) {
        sdWdx += ws_Sk_[j]*ws_dndx_[offSetDnDxIr+j];
      }

      // consistent mass for mu*I
      lhs[ir*nodesPerElement+ir] += rIr*extCoeffIp*ipFactor*rowSumLump_;

      // column ic;
      double rhsGalSum = 0.0;
      for ( int ic = 0; ic < nodesPerElement; ++ic) {
        // offset for column dndx
        const int offSetDnDxIc = nDim_*nodesPerElement*ip + ic*nDim_;
        const double rIc = ws_shape_function_[ipNpe+ic];

        double lhsSupgSum = 0.0;
        for ( int j = 0; j < nDim_; ++j ) {
          lhsSupgSum += ws_Sk_[j]*ws_dndx_[offSetDnDxIc+j];
        }

        // Galerkin
        lhs[ir*nodesPerElement+ic] += (-sdWdx + rIr*extCoeffIp*rIc*consistentMass_)*ipFactor;
        rhsGalSum += (-sdWdx + rIr*extCoeffIp*rIc*consistentMass_)*ws_intensity_[ic];

        // SUPG
        lhs[ir*nodesPerElement+ic] += tau*sdWdx*(lhsSupgSum + extCoeffIp*rIc)*ipFactor;
      }
      // Galerkin; (-sj I dwdxj + w extCoeff I_ip - w S_ip)*weight*detJ
      rhs[ir] -= (rhsGalSum - rIr*(radSrcIp + isotropicScatterIp - extCoeffIp*ws_intensity_[ir]*rowSumLump_))*ipFactor;

      // SUPG right hand side is: tau*sj*dWdxj*(Res)*weight*detJ
      rhs[ir] -= tau*sdWdx*residual*ipFactor;
    }
  }
}
  
} // namespace nalu
} // namespace Sierra
