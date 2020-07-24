/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleMomentumEdgeSolverAlgorithm.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#define _MODIFY_VELOCITY_ 0
#define _MODIFY_VISCOSITY_ 0
#define _AXISYMMETRIC_JET_ 1
#define _BACKWARD_FACING_STEP_ 0

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleMomentumEdgeSolverAlgorithm - add LHS/RHS for momentum
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeSolverAlgorithm::AssembleMomentumEdgeSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.does_mesh_move()),
    includeDivU_(realm_.get_divU()),
    velocityRTM_(NULL),
    velocity_(NULL),
    coordinates_(NULL),
    dudx_(NULL),
    density_(NULL),
    viscosity_(NULL),
    edgeAreaVec_(NULL),
    massFlowRate_(NULL),
    dualNodalVolume_(NULL),
    pecletFunction_(NULL),
    Cw_(realm.get_turb_model_constant(TM_Cw))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  // extract viscosity  name
  const std::string viscName = realm_.is_turbulent()
    ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  viscosity_lam_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  edgeAreaVec_ = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
  massFlowRate_ = meta_data.get_field<ScalarFieldType>(stk::topology::EDGE_RANK, "mass_flow_rate");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function<double>(velocity_->name());

  // sgs stress perturbation
  if(realm_.solutionOptions_->momentumPerturb_) {

    NaluEnv::self().naluOutputP0() << std::endl << "LES sgs stress perturbation model active:" << std::endl;

    // sgs stress magnitude perturbation
    if(realm_.solutionOptions_->momentumMagnitudePerturb_) {

      coeff_kk_ = realm_.solutionOptions_->momentumMagnitudePerturbCoeff_;

      NaluEnv::self().naluOutputP0() << "LES sgs magnitude perturbation: coeff_kk = " << coeff_kk_ << std::endl;

    }

    // sgs stress eigenvalue perturbation
    if(realm_.solutionOptions_->momentumEigenvaluePerturb_) {

      const int biasTowards = realm_.solutionOptions_->momentumEigenvaluePerturbBiasTowards_;

      if ( biasTowards == 1 ) {
        BinvXt_[0] =  2.0/3.0;
        BinvXt_[1] = -1.0/3.0;
        BinvXt_[2] = -1.0/3.0;
      }
      else if ( biasTowards == 2 ) {
        BinvXt_[0] =  1.0/6.0;
        BinvXt_[1] =  1.0/6.0;
        BinvXt_[2] = -1.0/3.0;
      }
      else if ( biasTowards == 3 ) {
        BinvXt_[0] = 0.0;
        BinvXt_[1] = 0.0;
        BinvXt_[2] = 0.0;
      }
      else {
        throw std::runtime_error("LES sgs stress eigenvalue perturbation: target eigenvalue incorrect");
      }

      deltaB_ = realm_.solutionOptions_->momentumEigenvaluePerturbDelta_;

      NaluEnv::self().naluOutputP0() << "LES sgs stress eigenvalue perturbation: towards = " << biasTowards << ", deltaB = " << deltaB_ << std::endl;

    }

    // sgs stress eigenvector perturbation
    if(realm_.solutionOptions_->momentumEigenvectorPerturb_) {

      eigenvectorPermutation_ = realm_.solutionOptions_->momentumEigenvectorPerturbPermutation_;

      if ( ! ( (eigenvectorPermutation_ == 1) || (eigenvectorPermutation_ == 2) || (eigenvectorPermutation_ == 3) ) ) {

        throw std::runtime_error("LES sgs stress eigenvector perturbation: eigenvector permutation incorrect");

      }

      NaluEnv::self().naluOutputP0() << "LES sgs stress eigenvector perturbation: permutation = " << eigenvectorPermutation_ << std::endl;

    }

    NaluEnv::self().naluOutputP0() << std::endl;

  }

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeSolverAlgorithm::~AssembleMomentumEdgeSolverAlgorithm()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const double invNdim = 1.0/double(nDim);

  const double small = 1.0e-16;

  // extract user advection options (allow to potentially change over time)
  const std::string dofName = "velocity";
  const double alpha = realm_.get_alpha_factor(dofName);
  const double alphaUpw = realm_.get_alpha_upw_factor(dofName);
  const double hoUpwind = realm_.get_upw_factor(dofName);
  const bool useLimiter = realm_.primitive_uses_limiter(dofName);

  // one minus flavor
  const double om_alpha = 1.0-alpha;
  const double om_alphaUpw = 1.0-alphaUpw;

  // space for LHS/RHS; always edge connectivity
  const int nodesPerEdge = 2;
  const int lhsSize = nDim*nodesPerEdge*nDim*nodesPerEdge;
  const int rhsSize = nDim*nodesPerEdge;
  std::vector<double> lhs(lhsSize);
  std::vector<double> rhs(rhsSize);
  std::vector<int> scratchIds(rhsSize);
  std::vector<double> scratchVals(rhsSize);
  std::vector<stk::mesh::Entity> connected_nodes(2);

  // area vector; gather into
  std::vector<double> areaVec(nDim);

  // pointer for fast access
  double *p_lhs = &lhs[0];
  double *p_rhs = &rhs[0];
  double *p_areaVec = &areaVec[0];

  // space for dui/dxj. This variable is the modified gradient with NOC
  std::vector<double> duidxj(nDim*nDim);

  // extrapolated value from the L/R direction 
  std::vector<double> uIpL(nDim);
  std::vector<double> uIpR(nDim);
  // limiter values from the L/R direction, 0:1
  std::vector<double> limitL(nDim,1.0); 
  std::vector<double> limitR(nDim,1.0);
  // extrapolated gradient from L/R direction
  std::vector<double> duL(nDim);
  std::vector<double> duR(nDim);
  
  // pointers for fast access
  double *p_duidxj = &duidxj[0];
  double *p_uIpL = &uIpL[0];
  double *p_uIpR = &uIpR[0];
  double *p_limitL = &limitL[0];
  double *p_limitR = &limitR[0];
  double *p_duL = &duL[0];
  double *p_duR = &duR[0];

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& edge_buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to edge area vector and mdot
    const double * av = stk::mesh::field_data(*edgeAreaVec_, b);
    const double * mdot = stk::mesh::field_data(*massFlowRate_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zeroing of lhs/rhs
      for ( int i = 0; i < lhsSize; ++i ) {
        p_lhs[i] = 0.0;
      }
      for ( int i = 0; i < rhsSize; ++i ) {
        p_rhs[i] = 0.0;
      }

      stk::mesh::Entity const * edge_node_rels = b.begin_nodes(k);

      // pointer to edge area vector
      for ( int j = 0; j < nDim; ++j )
        p_areaVec[j] = av[k*nDim+j];
      const double tmdot = mdot[k];

      // sanity check on number or nodes
      ThrowAssert( b.num_nodes(k) == 2 );

      // left and right nodes
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];

      connected_nodes[0] = nodeL;
      connected_nodes[1] = nodeR;

      // extract nodal fields
      const double * coordL = stk::mesh::field_data(*coordinates_, nodeL);
      const double * coordR = stk::mesh::field_data(*coordinates_, nodeR);

      const double * dudxL = stk::mesh::field_data(*dudx_, nodeL);
      const double * dudxR = stk::mesh::field_data(*dudx_, nodeR);

      const double * vrtmL = stk::mesh::field_data(*velocityRTM_, nodeL);
      const double * vrtmR = stk::mesh::field_data(*velocityRTM_, nodeR);

#if _MODIFY_VELOCITY_
      double * uNp1L = stk::mesh::field_data(velocityNp1, nodeL);
      double * uNp1R = stk::mesh::field_data(velocityNp1, nodeR);

      // Modify velocity at the top boundary to enforce outflow

#if _AXISYMMETRIC_JET_
      //if( 0.5*( coordR[2] + coordL[2] ) > 0.75 ) {
      //if( 0.5*( coordR[2] + coordL[2] ) > 0.7 ) {
      if( 0.5*( coordR[2] + coordL[2] ) > 0.675 ) {

        const double vel_x = 0.0;
        const double vel_y = 0.0;
        const double vel_z = 0.92775;

        //const double original_weight_factor = 0.999;
        const double original_weight_factor = 0.995;


        uNp1L[0] = original_weight_factor*uNp1L[0] + (1.0 - original_weight_factor)*vel_x;
        uNp1R[0] = original_weight_factor*uNp1R[0] + (1.0 - original_weight_factor)*vel_x;

        uNp1L[1] = original_weight_factor*uNp1L[1] + (1.0 - original_weight_factor)*vel_y;
        uNp1R[1] = original_weight_factor*uNp1R[1] + (1.0 - original_weight_factor)*vel_y;

        uNp1L[2] = original_weight_factor*uNp1L[2] + (1.0 - original_weight_factor)*vel_z;
        uNp1R[2] = original_weight_factor*uNp1R[2] + (1.0 - original_weight_factor)*vel_z;

      }
#endif

#else
      const double * uNp1L = stk::mesh::field_data(velocityNp1, nodeL);
      const double * uNp1R = stk::mesh::field_data(velocityNp1, nodeR);
#endif

      const double densityL = *stk::mesh::field_data(densityNp1, nodeL);
      const double densityR = *stk::mesh::field_data(densityNp1, nodeR);

#if _MODIFY_VISCOSITY_
      double viscosityL = *stk::mesh::field_data(*viscosity_, nodeL);
      double viscosityR = *stk::mesh::field_data(*viscosity_, nodeR);

      double viscosityLamL = *stk::mesh::field_data(*viscosity_lam_, nodeL);
      double viscosityLamR = *stk::mesh::field_data(*viscosity_lam_, nodeR);

      // Modify molecular viscosity at the boundaries to laminarize flow

#if _AXISYMMETRIC_JET_
      //if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) > 0.20 ) || ( 0.5*( coordR[2] + coordL[2] ) > 0.575 ) ) {
      //if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) > 0.20 ) || ( 0.5*( coordR[2] + coordL[2] ) > 0.675 ) ) {
      //if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) > 0.23 ) || ( 0.5*( coordR[2] + coordL[2] ) > 0.75 ) ) {
      if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) > 0.23 ) || ( 0.5*( coordR[2] + coordL[2] ) > 0.70 ) ) {

        // Increase viscosity by 100x
        //double increased_viscosity = 10.0*1.8e-5;
        //double increased_viscosity = 20.0*1.8e-5;
        //double increased_viscosity = 100.0*1.8e-5;
        //double increased_viscosity = 1000.0*1.8e-5;

        double increased_viscosity = 1.0e-2*1.8e-5;

        viscosityLamL = increased_viscosity;
        viscosityLamR = increased_viscosity;

      }
#endif

#if _BACKWARD_FACING_STEP_
      if( 0.5*( coordR[0] + coordL[0] ) > 26.0 ) {

        // Increase viscosity by 100x
        const double increased_viscosity = 100.0*0.000196078;

        viscosityLamL = increased_viscosity;
        viscosityLamR = increased_viscosity;

      }
#endif

      // Switch turbulent viscosity off near the jet injection
#if 0
//#if _AXISYMMETRIC_JET_
      if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) < 0.0145 ) && ( 0.5*( coordR[2] + coordL[2] ) < 0.026 ) ) {
      //if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) < 0.0145 ) && ( 0.5*( coordR[2] + coordL[2] ) < 0.13 ) ) {
      //if( ( sqrt( pow( 0.5*( coordR[0] + coordL[0] ), 2.0 ) + pow( 0.5*( coordR[1] + coordL[1] ), 2.0 ) ) < 0.0145 ) && ( 0.5*( coordR[2] + coordL[2] ) < 0.195 ) ) {

        viscosityL = 1.8e-5;
        viscosityR = 1.8e-5;

      }
#endif

#else
      const double viscosityL = *stk::mesh::field_data(*viscosity_, nodeL);
      const double viscosityR = *stk::mesh::field_data(*viscosity_, nodeR);

      const double viscosityLamL = *stk::mesh::field_data(*viscosity_lam_, nodeL);
      const double viscosityLamR = *stk::mesh::field_data(*viscosity_lam_, nodeR);
#endif

      const double volL = *stk::mesh::field_data( *dualNodalVolume_, nodeL );
      const double volR = *stk::mesh::field_data( *dualNodalVolume_, nodeR );

      // copy in extrapolated values
      for ( int i = 0; i < nDim; ++i ) {
        // extrapolated du
        p_duL[i] = 0.0;
        p_duR[i] = 0.0;
        const int offSet = nDim*i;
        for ( int j = 0; j < nDim; ++j ) {
          const double dxj = 0.5*(coordR[j] - coordL[j]);
          p_duL[i] += dxj*dudxL[offSet+j];
          p_duR[i] += dxj*dudxR[offSet+j];
        }
      }

      // compute geometry
      double axdx = 0.0;
      double asq = 0.0;
      double udotx = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = p_areaVec[j];
        const double dxj = coordR[j] - coordL[j];
        axdx += axj*dxj;
        asq += axj*axj;
        udotx += 0.5*dxj*(vrtmL[j] + vrtmR[j]);
      }

      const double inv_axdx = 1.0/axdx;

      // ip props
      const double rhoIp = 0.5*(densityL + densityR);
      const double viscIp = 0.5*(viscosityL + viscosityR);
      const double viscLamIp = 0.5*(viscosityLamL + viscosityLamR);
      const double diffIp = 0.5*(viscosityL/densityL + viscosityR/densityR);
      const double volIp = 0.5*(volL + volR);

      // Peclet factor
      const double pecfac = pecletFunction_->execute(std::abs(udotx)/(diffIp+small));
      const double om_pecfac = 1.0-pecfac;

      // determine limiter if applicable
      if ( useLimiter ) {
        for ( int i = 0; i < nDim; ++i ) {
          const double dq = uNp1R[i] - uNp1L[i];
          const double dqMl = 2.0*2.0*p_duL[i] - dq;
          const double dqMr = 2.0*2.0*p_duR[i] - dq;
          p_limitL[i] = van_leer(dqMl, dq, small);
          p_limitR[i] = van_leer(dqMr, dq, small);
        }
      }

      // final upwind extrapolation; with limiter
      for ( int i = 0; i < nDim; ++i ) {
        p_uIpL[i] = uNp1L[i] + p_duL[i]*hoUpwind*p_limitL[i];
        p_uIpR[i] = uNp1R[i] - p_duR[i]*hoUpwind*p_limitR[i];
      }

      /*
        form duidxj with over-relaxed procedure of Jasak:

        dui/dxj = GjUi +[(uiR - uiL) - GlUi*dxl]*Aj/AxDx
        where Gp is the interpolated pth nodal gradient for ui
      */
      for ( int i = 0; i < nDim; ++i ) {

        // difference between R and L nodes for component i
        const double uidiff = uNp1R[i] - uNp1L[i];

        // offset into all forms of dudx
        const int offSetI = nDim*i;

        // start sum for NOC contribution
        double GlUidxl = 0.0;
        for ( int l = 0; l< nDim; ++l ) {
          const int offSetIL = offSetI+l;
          const double dxl = coordR[l] - coordL[l];
          const double GlUi = 0.5*(dudxL[offSetIL] + dudxR[offSetIL]);
          GlUidxl += GlUi*dxl;
        }

        // form full tensor dui/dxj with NOC
        for ( int j = 0; j < nDim; ++j ) {
          const int offSetIJ = offSetI+j;
          const double axj = p_areaVec[j];
          const double GjUi = 0.5*(dudxL[offSetIJ] + dudxR[offSetIJ]);
          p_duidxj[offSetIJ] = GjUi + (uidiff - GlUidxl)*axj*inv_axdx;
        }
      }

      // sgs stress perturbation
      if(realm_.solutionOptions_->momentumPerturb_) {

        // calculate resolved_kk
        double uiIp_[3];
        uiIp_[0] = 0.5*(uNp1R[0] + uNp1L[0]);
        uiIp_[1] = 0.5*(uNp1R[1] + uNp1L[1]);
        uiIp_[2] = 0.5*(uNp1R[2] + uNp1L[2]);

        resolved_kk = uiIp_[0]*uiIp_[0] + uiIp_[1]*uiIp_[1] + uiIp_[2]*uiIp_[2] + 1.0e-20;
        //NaluEnv::self().naluOutputP0() << "resolved_kk: " << resolved_kk << std::endl;

        // compute the resolved strain tensor S_
        for ( int i = 0; i < nDim; ++i ) {
          const int offSetI = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            S_[i][j] = 0.5*(p_duidxj[offSetI+j] + p_duidxj[offSetTrans]);
          }
        }
        //NaluEnv::self().naluOutputP0() << "Trace S_ij: " << S_[0][0] + S_[1][1] + S_[2][2] << std::endl;

        // correct non-zero trace S_ tensor. Even for incompressible flow simulation, the edge-based velocity field may not be discretely divergence-free. The eigenvalues are slightly modified but not the eigenvectors.
        const double tmp = (S_[0][0] + S_[1][1] + S_[2][2])/3.0;
        S_[0][0] -= tmp;
        S_[1][1] -= tmp;
        S_[2][2] -= tmp;
        //NaluEnv::self().naluOutputP0() << "Trace S_ij: " << S_[0][0] + S_[1][1] + S_[2][2] << std::endl;

        // compute filter size
        filter = std::pow(volIp, invNdim);
        //NaluEnv::self().naluOutputP0() << "Filter size: " << filter << std::endl;

        // compute nu_sgs
        nu_sgs = (viscIp - viscLamIp)/rhoIp;
        //NaluEnv::self().naluOutputP0() << "nu: " << viscIp/rhoIp << "   nu_lam: " << viscLamIp/rhoIp << "   nu_sgs: " << nu_sgs << std::endl;

#if 0
        // compute S_mag
        const double S_mag = sqrt( 2.0*( S_[0][0]*S_[0][0] +
                                         S_[1][1]*S_[1][1] +
                                         S_[2][2]*S_[2][2] ) +
                                   4.0*( S_[0][1]*S_[0][1] +
                                         S_[1][2]*S_[1][2] +
                                         S_[2][0]*S_[2][0] ) );
        //NaluEnv::self().naluOutputP0() << "S_mag: " << S_mag << std::endl;

        // calculate modeled_kk (A. Yoshizawa. Phys. Fluids 29, 2152, 1986; Moin et al. Phys. Fluids A, 3 (11), 1991)
        //const double C_I = Cw_*Cw_; // using the same coefficient as the sgs WALE model
        const double C_I = 0.5*sqrt(3.0)*Cw_*Cw_; // using the same coefficient as the sgs WALE model
        modeled_kk = 2.0*C_I*filter*filter*S_mag*S_mag + 1.0e-20;
#else
        // compute S_mag
        const double S_mag = sqrt( 0.5*( 4.0*S_[0][0]*S_[0][0] + 4.0*S_[0][1]*S_[0][1] + 4.0*S_[0][2]*S_[0][2]
                                       + 4.0*S_[1][0]*S_[1][0] + 4.0*S_[1][1]*S_[1][1] + 4.0*S_[1][2]*S_[1][2]
                                       + 4.0*S_[2][0]*S_[2][0] + 4.0*S_[2][1]*S_[2][1] + 4.0*S_[2][2]*S_[2][2] ) );
        //NaluEnv::self().naluOutputP0() << "S_mag: " << S_mag << std::endl;

        // calculate modeled_kk
        modeled_kk = sqrt(3.0)*nu_sgs*S_mag + 1.0e-20;
        //modeled_kk = sqrt(3.0)*(0.5*0.5*nu_sgs/(Cw_*Cw_))*S_mag + 1.0e-20; // Modification to recover Nicoud et al. original Cw_ value (Cw_ = 0.5).
        //modeled_kk = 2.0*0.0886*filter*filter*S_mag*S_mag + 1.0e-20; // Yoshizawa model proposed in Vreman et al. 1994.
#endif
        //NaluEnv::self().naluOutputP0() << "modeled_kk: " << modeled_kk << std::endl;

        // calculate total_kk = resolved_kk + modeled_kk
        total_kk = resolved_kk + modeled_kk + 1.0e-20;
        //NaluEnv::self().naluOutputP0() << "total_kk: " << total_kk << std::endl;

        // compute the modeled anisotropy sgs stress tensor A_
        for ( int i = 0; i < nDim; ++i ) {
          for ( int j = 0; j < nDim; ++j ) {
            A_[i][j] = (-2.0)*(nu_sgs/total_kk)*S_[i][j];
          }
        }
        //if( A_[0][0] < (-1.0/3.0) || A_[0][0] > (2.0/3.0) || A_[1][1] < (-1.0/3.0) || A_[1][1] > (2.0/3.0) || A_[2][2] < (-1.0/3.0) || A_[2][2] > (2.0/3.0) || A_[0][1] < -0.5 || A_[0][1] > 0.5 || A_[0][2] < -0.5 || A_[0][2] > 0.5 || A_[1][0] < -0.5 || A_[1][0] > 0.5 || A_[1][2] < -0.5 || A_[1][2] > 0.5 || A_[2][0] < -0.5 || A_[2][0] > 0.5 || A_[2][1] < -0.5 || A_[2][1] > 0.5 ) { NaluEnv::self().naluOutputP0() << "a11: " << A_[0][0] << "  a12: " << A_[0][1] << "  a13: " << A_[0][2] << "  a21: " << A_[1][0] << "  a22: " << A_[1][1] << "  a23: " << A_[1][2] << "  a31: " << A_[2][0] << "  a32: " << A_[2][1] << "  a33: " << A_[2][2] << "   total_kk: " << total_kk << std::endl; } else { NaluEnv::self().naluOutputP0() << "total_kk: " << total_kk << std::endl;} 

        // perform the decomposition
        diagonalize(A_, Q_, D_);
        //NaluEnv::self().naluOutputP0() << "Lambda1: " << D_[0][0] << "   Lambda2: " << D_[1][1] << "   Lambda3: " << D_[2][2] << "   Sum: " << D_[0][0] + D_[1][1] + D_[2][2] << "   Sum - A_trace: " << (D_[0][0] + D_[1][1] + D_[2][2]) - (A_[0][0] + A_[1][1] + A_[2][2]) << std::endl;
        //NaluEnv::self().naluOutputP0() << D_[0][0] << "  " << D_[0][1] << "  " << D_[0][2] << "  " << D_[1][0] << "  " << D_[1][1] << "  " << D_[1][2] << "  " << D_[2][0] << "  " <<D_[2][1] << "  " << D_[2][2] << std::endl;

        // sort Q and D: eigenvectors and eigenvalues need to be reorganized at the same time to recover the same tensor
        sort(Q_, D_);
        //NaluEnv::self().naluOutputP0() << "Lambda1: " << D_[0][0] << "   Lambda2: " << D_[1][1] << "   Lambda3: " << D_[2][2] << "   Sum: " << D_[0][0] + D_[1][1] + D_[2][2] << std::endl;
        //NaluEnv::self().naluOutputP0() <<"total_kk: " << total_kk << "  nu_sgs: " << nu_sgs << "  x0: " << 0.0*(D_[0][0] - D_[1][1]) + 1.0*(2.0*D_[1][1] - 2.0*D_[2][2]) + 0.5*(3.0*D_[2][2] + 1.0) << "  x1: " << 0.0*(D_[0][0] - D_[1][1]) + 0.0*(2.0*D_[1][1] - 2.0*D_[2][2]) + (sqrt(3.0)/2.0)*(3.0*D_[2][2] + 1.0) << std::endl;

        // perturb
        perturb(Q_,D_);
        //NaluEnv::self().naluOutputP0() << "Lambda1: " << D_[0][0] << "   Lambda2: " << D_[1][1] << "   Lambda3: " << D_[2][2] << "   Sum: " << D_[0][0] + D_[1][1] + D_[2][2] << std::endl;

        // form perturbed stress tensor (deviatoric term absorbed into pressure)
        form_perturbed_stress(D_, Q_, A_);

        // sgs stress magnitude perturbation
        if(realm_.solutionOptions_->momentumMagnitudePerturb_) {

          // -resolved_kk <= modeled_kk + delta_kk <= total_kk
          // Therefore, -resolved_kk - modeled_kk <= delta_kk <= total_kk - modeled_kk

          const double delta_kk_ = (coeff_kk_ > 0.0) ? (coeff_kk_*resolved_kk) : (coeff_kk_*(resolved_kk + modeled_kk));
          //NaluEnv::self().naluOutputP0() << "delta_kk: " << delta_kk << std::endl;

          modeled_kk += delta_kk_;
          //NaluEnv::self().naluOutputP0() << "modeled_kk: " << modeled_kk << std::endl;

          // calculate total_kk = resolved_kk + modeled_kk
          total_kk = resolved_kk + modeled_kk + 1.0e-20;
          //NaluEnv::self().naluOutputP0() << "total_kk: " << total_kk << std::endl;

        }

        // In the case of eigenvector perturbation with permutation 3 -> ensure that viscous dissipation is larger than SGS dissipation
        if(realm_.solutionOptions_->momentumEigenvectorPerturb_) {

          if(eigenvectorPermutation_ == 3) {

            // SGS and viscous dissipations
            double sgsDissipation = 0.0;
            double visDissipation = 0.0;
            for ( int i = 0; i < nDim; ++i ) {
              for ( int j = 0; j < nDim; ++j ) {
                sgsDissipation += (total_kk*A_[i][j])*S_[i][j];
                visDissipation += 2.0*(viscLamIp/rhoIp)*S_[i][j]*S_[i][j];
              }
            }

            // scaling factor
            const double scalingFactor = fabs(visDissipation/sgsDissipation);
            //NaluEnv::self().naluOutputP0() << "Scaling factor: " << scalingFactor << std::endl;

            // scale sgs (total) production if it is larger than viscous dissipation
            if(scalingFactor < 1.0) {

              //total_kk *= scalingFactor;
              total_kk *= scalingFactor;
              total_kk += 1.0e-20;
              //NaluEnv::self().naluOutputP0() << "total_kk: " << total_kk << std::endl;

            }

          }

        }

#if 0 
        // Check if the resolved + perturbed_modeled tensor is realizable

        double anisotropyNonFilterAdvectionTensor_[3][3];

        // compute the anisotropy nonlinear filtered advection term
        for ( int i = 0; i < nDim; ++i ) {
          const int offSetI = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            anisotropyNonFilterAdvectionTensor_[i][j] = (1.0/total_kk)*( uiIp_[i]*uiIp_[j] + total_kk*A_[i][j] - ( (i ==j) ? (resolved_kk/3.0) : 0.0 ) );
          }
        }

        // correct non-zero trace anisotropyNonFilterAdvectionTensor_
        const double tmpAuxAux = (anisotropyNonFilterAdvectionTensor_[0][0] + anisotropyNonFilterAdvectionTensor_[1][1] + anisotropyNonFilterAdvectionTensor_[2][2])/3.0;
        anisotropyNonFilterAdvectionTensor_[0][0] -= tmpAuxAux;
        anisotropyNonFilterAdvectionTensor_[1][1] -= tmpAuxAux;
        anisotropyNonFilterAdvectionTensor_[2][2] -= tmpAuxAux;

        //if( anisotropyNonFilterAdvectionTensor_[0][0] < (-1.0/3.0) || anisotropyNonFilterAdvectionTensor_[0][0] > (2.0/3.0) || anisotropyNonFilterAdvectionTensor_[1][1] < (-1.0/3.0) || anisotropyNonFilterAdvectionTensor_[1][1] > (2.0/3.0) || anisotropyNonFilterAdvectionTensor_[2][2] < (-1.0/3.0) || anisotropyNonFilterAdvectionTensor_[2][2] > (2.0/3.0) || anisotropyNonFilterAdvectionTensor_[0][1] < -0.5 || anisotropyNonFilterAdvectionTensor_[0][1] > 0.5 || anisotropyNonFilterAdvectionTensor_[0][2] < -0.5 || anisotropyNonFilterAdvectionTensor_[0][2] > 0.5 || anisotropyNonFilterAdvectionTensor_[1][0] < -0.5 || anisotropyNonFilterAdvectionTensor_[1][0] > 0.5 || anisotropyNonFilterAdvectionTensor_[1][2] < -0.5 || anisotropyNonFilterAdvectionTensor_[1][2] > 0.5 || anisotropyNonFilterAdvectionTensor_[2][0] < -0.5 || anisotropyNonFilterAdvectionTensor_[2][0] > 0.5 || anisotropyNonFilterAdvectionTensor_[2][1] < -0.5 || anisotropyNonFilterAdvectionTensor_[2][1] > 0.5 ) { NaluEnv::self().naluOutputP0() << "a11: " << anisotropyNonFilterAdvectionTensor_[0][0] << "  a12: " << anisotropyNonFilterAdvectionTensor_[0][1] << "  a13: " << anisotropyNonFilterAdvectionTensor_[0][2] << "  a21: " << anisotropyNonFilterAdvectionTensor_[1][0] << "  a22: " << anisotropyNonFilterAdvectionTensor_[1][1] << "  a23: " << anisotropyNonFilterAdvectionTensor_[1][2] << "  a31: " << anisotropyNonFilterAdvectionTensor_[2][0] << "  a32: " << anisotropyNonFilterAdvectionTensor_[2][1] << "  a33: " << anisotropyNonFilterAdvectionTensor_[2][2] << std::endl; }

#endif

      }

      // lhs diffusion; only -mu*dui/dxj*Aj contribution for now
      const double dlhsfac = -viscIp*asq*inv_axdx;

      for ( int i = 0; i < nDim; ++i ) {

        // 2nd order central
        const double uiIp = 0.5*(uNp1R[i] + uNp1L[i]);

        // upwind
        const double uiUpwind = (tmdot > 0) ? alphaUpw*p_uIpL[i] + om_alphaUpw*uiIp
          : alphaUpw*p_uIpR[i] + om_alphaUpw*uiIp;

        // generalized central (2nd and 4th order)
        const double uiHatL = alpha*p_uIpL[i] + om_alpha*uiIp;
        const double uiHatR = alpha*p_uIpR[i] + om_alpha*uiIp;
        const double uiCds = 0.5*(uiHatL + uiHatR);

        // total advection; pressure contribution in time term expression
        const double aflux = tmdot*(pecfac*uiUpwind + om_pecfac*uiCds);

        // divU
        double divU = 0.0;
        for ( int j = 0; j < nDim; ++j)
          divU += p_duidxj[j*nDim+j];

        // diffusive flux; viscous tensor doted with area vector
        double dflux = 2.0/3.0*viscIp*divU*p_areaVec[i]*includeDivU_;
        const int offSetI = nDim*i;
        if(realm_.solutionOptions_->momentumPerturb_) {
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            const double axj = p_areaVec[j];

            // diffusion flux
            dflux += -viscLamIp*(p_duidxj[offSetI+j] + p_duidxj[offSetTrans])*axj;

            // sgs stress tensor flux
            dflux += rhoIp*total_kk*A_[i][j]*axj;
          }
        } else {
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            const double axj = p_areaVec[j];
            dflux += -viscIp*(p_duidxj[offSetI+j] + p_duidxj[offSetTrans])*axj;
          }
        }

        // residual for total flux
        const double tflux = aflux + dflux;
        const int indexL = i;
        const int indexR = i + nDim;

        // total flux left
        p_rhs[indexL] -= tflux;
        // total flux right
        p_rhs[indexR] += tflux;

        // setup for LHS
        const int rowL = indexL * nodesPerEdge*nDim;
        const int rowR = indexR * nodesPerEdge*nDim;

        //==============================
        // advection first
        //==============================
        const int rLiL = rowL+indexL;
        const int rLiR = rowL+indexR;
        const int rRiL = rowR+indexL;
        const int rRiR = rowR+indexR;

        // upwind advection (includes 4th); left node
        double alhsfac = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[rLiL] += alhsfac;
        p_lhs[rRiL] -= alhsfac;

        // upwind advection (incldues 4th); right node
        alhsfac = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[rRiR] -= alhsfac;
        p_lhs[rLiR] += alhsfac;

        // central; left; collect terms on alpha and alphaUpw
        alhsfac = 0.5*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);
        p_lhs[rLiL] += alhsfac;
        p_lhs[rLiR] += alhsfac;
        // central; right
        p_lhs[rRiL] -= alhsfac;
        p_lhs[rRiR] -= alhsfac;

        //==============================
        // diffusion second
        //==============================
        const double axi = p_areaVec[i];

        //diffusion; row IL
        p_lhs[rLiL] -= dlhsfac;
        p_lhs[rLiR] += dlhsfac;

        // diffusion; row IR
        p_lhs[rRiL] += dlhsfac;
        p_lhs[rRiR] -= dlhsfac;

        // more diffusion; see theory manual
        for ( int j = 0; j < nDim; ++j ) {
          const double lhsfacNS = -viscIp*axi*p_areaVec[j]*inv_axdx;

          const int colL = j;
          const int colR = j + nDim;

          // first left; IL,IL; IL,IR
          p_lhs[rowL + colL] -= lhsfacNS;
          p_lhs[rowL + colR] += lhsfacNS;

          // now right, IR,IL; IR,IR
          p_lhs[rowR + colL] += lhsfacNS;
          p_lhs[rowR + colR] -= lhsfacNS;
        }

      }
      
      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }
  }
}

//--------------------------------------------------------------------------
//-------- diagonalize -----------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::diagonalize(
  const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3])
{
  /*  
    obtained from: 
    http://stackoverflow.com/questions/4372224/
    fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition

    A must be a symmetric matrix.
    returns Q and D such that 
    Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
  */

  const int maxsteps=24;
  int k0, k1, k2;
  double o[3], m[3];
  double q [4] = {0.0,0.0,0.0,1.0};
  double jr[4];
  double sqw, sqx, sqy, sqz;
  double tmp1, tmp2, mq;
  double AQ[3][3];
  double thet, sgn, t, c;
  for(int i=0;i < maxsteps;++i) {
    // quat to matrix
    sqx      = q[0]*q[0];
    sqy      = q[1]*q[1];
    sqz      = q[2]*q[2];
    sqw      = q[3]*q[3];
    Q[0][0]  = ( sqx - sqy - sqz + sqw);
    Q[1][1]  = (-sqx + sqy - sqz + sqw);
    Q[2][2]  = (-sqx - sqy + sqz + sqw);
    tmp1     = q[0]*q[1];
    tmp2     = q[2]*q[3];
    Q[1][0]  = 2.0 * (tmp1 + tmp2);
    Q[0][1]  = 2.0 * (tmp1 - tmp2);
    tmp1     = q[0]*q[2];
    tmp2     = q[1]*q[3];
    Q[2][0]  = 2.0 * (tmp1 - tmp2);
    Q[0][2]  = 2.0 * (tmp1 + tmp2);
    tmp1     = q[1]*q[2];
    tmp2     = q[0]*q[3];
    Q[2][1]  = 2.0 * (tmp1 + tmp2);
    Q[1][2]  = 2.0 * (tmp1 - tmp2);

    // AQ = A * Q
    AQ[0][0] = Q[0][0]*A[0][0]+Q[1][0]*A[0][1]+Q[2][0]*A[0][2];
    AQ[0][1] = Q[0][1]*A[0][0]+Q[1][1]*A[0][1]+Q[2][1]*A[0][2];
    AQ[0][2] = Q[0][2]*A[0][0]+Q[1][2]*A[0][1]+Q[2][2]*A[0][2];
    AQ[1][0] = Q[0][0]*A[0][1]+Q[1][0]*A[1][1]+Q[2][0]*A[1][2];
    AQ[1][1] = Q[0][1]*A[0][1]+Q[1][1]*A[1][1]+Q[2][1]*A[1][2];
    AQ[1][2] = Q[0][2]*A[0][1]+Q[1][2]*A[1][1]+Q[2][2]*A[1][2];
    AQ[2][0] = Q[0][0]*A[0][2]+Q[1][0]*A[1][2]+Q[2][0]*A[2][2];
    AQ[2][1] = Q[0][1]*A[0][2]+Q[1][1]*A[1][2]+Q[2][1]*A[2][2];
    AQ[2][2] = Q[0][2]*A[0][2]+Q[1][2]*A[1][2]+Q[2][2]*A[2][2];
    // D = Qt * AQ
    D[0][0] = AQ[0][0]*Q[0][0]+AQ[1][0]*Q[1][0]+AQ[2][0]*Q[2][0];
    D[0][1] = AQ[0][0]*Q[0][1]+AQ[1][0]*Q[1][1]+AQ[2][0]*Q[2][1];
    D[0][2] = AQ[0][0]*Q[0][2]+AQ[1][0]*Q[1][2]+AQ[2][0]*Q[2][2];
    D[1][0] = AQ[0][1]*Q[0][0]+AQ[1][1]*Q[1][0]+AQ[2][1]*Q[2][0];
    D[1][1] = AQ[0][1]*Q[0][1]+AQ[1][1]*Q[1][1]+AQ[2][1]*Q[2][1];
    D[1][2] = AQ[0][1]*Q[0][2]+AQ[1][1]*Q[1][2]+AQ[2][1]*Q[2][2];
    D[2][0] = AQ[0][2]*Q[0][0]+AQ[1][2]*Q[1][0]+AQ[2][2]*Q[2][0];
    D[2][1] = AQ[0][2]*Q[0][1]+AQ[1][2]*Q[1][1]+AQ[2][2]*Q[2][1];
    D[2][2] = AQ[0][2]*Q[0][2]+AQ[1][2]*Q[1][2]+AQ[2][2]*Q[2][2];
    o[0]    = D[1][2];
    o[1]    = D[0][2];
    o[2]    = D[0][1];
    m[0]    = std::abs(o[0]);
    m[1]    = std::abs(o[1]);
    m[2]    = std::abs(o[2]);

    k0      = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
    k1      = (k0+1)%3;
    k2      = (k0+2)%3;
    if (o[k0]==0.0) {
      break;  // diagonal already
    }
    thet    = (D[k2][k2]-D[k1][k1])/(2.0*o[k0]);
    sgn     = (thet > 0.0)?1.0:-1.0;
    thet   *= sgn; // make it positive
    t       = sgn /(thet +((thet < 1.E6)? std::sqrt(thet*thet+1.0):thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
    c       = 1.0/std::sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c 
    if(c==1.0) {
      break;  // no room for improvement - reached machine precision.
    }
    jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
    jr[k0]  = sgn*std::sqrt((1.0-c)/2.0);  // using 1/2 angle identity sin(a/2) = std::sqrt((1-cos(a))/2)  
    jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
    jr[3 ]  = std::sqrt(1.0f - jr[k0] * jr[k0]);
    if(jr[3]==1.0) {
      break; // reached limits of floating point precision
    }
    q[0]    = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
    q[1]    = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
    q[2]    = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
    q[3]    = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
    mq      = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0]   /= mq;
    q[1]   /= mq;
    q[2]   /= mq;
    q[3]   /= mq;
  }
}

//--------------------------------------------------------------------------
//-------- sort ------------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::sort(
  double (&Q)[3][3], double (&D)[3][3])
{
  // Goal: sort diagonalization eigenvalues from high to low; save off row in D from high to low; reorder Q tensor accordingly

  for(int i = 0; i < 2; ++i) {

    if(D[0][0] < D[1][1]) {

      double tempEigenvalue = D[0][0];
      D[0][0] = D[1][1];
      D[1][1] = tempEigenvalue;

      double tempEigenvector0 = Q[0][0]; double tempEigenvector1 = Q[1][0]; double tempEigenvector2 = Q[2][0];
      Q[0][0] = Q[0][1]; Q[1][0] = Q[1][1]; Q[2][0] = Q[2][1];
      Q[0][1] = tempEigenvector0; Q[1][1] = tempEigenvector1; Q[2][1] = tempEigenvector2;

    }

    if(D[1][1] < D[2][2]) {

      double tempEigenvalue = D[1][1];
      D[1][1] = D[2][2];
      D[2][2] = tempEigenvalue;

      double tempEigenvector0 = Q[0][1]; double tempEigenvector1 = Q[1][1]; double tempEigenvector2 = Q[2][1];
      Q[0][1] = Q[0][2]; Q[1][1] = Q[1][2]; Q[2][1] = Q[2][2];
      Q[0][2] = tempEigenvector0; Q[1][2] = tempEigenvector1; Q[2][2] = tempEigenvector2;

    }

  }
}

//--------------------------------------------------------------------------
//-------- perturb ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::perturb(
  double (&Q)[3][3], double (&D)[3][3])
{
  // sgs stress eigenvalue perturbation
  if(realm_.solutionOptions_->momentumEigenvaluePerturb_) {

    // extract sorted by size (L1 > L2 > L3)
    const double Lambda1 = D[0][0];
    const double Lambda2 = D[1][1];
    const double Lambda3 = D[2][2];

    const double pLambda1 = (1.0 - deltaB_)*Lambda1 + deltaB_*BinvXt_[0];
    const double pLambda2 = (1.0 - deltaB_)*Lambda2 + deltaB_*BinvXt_[1];
    const double pLambda3 = (1.0 - deltaB_)*Lambda3 + deltaB_*BinvXt_[2];

    D[0][0] = pLambda1;
    D[1][1] = pLambda2;
    D[2][2] = pLambda3;

  }

  // sgs stress eigenvector perturbation
  if(realm_.solutionOptions_->momentumEigenvectorPerturb_) {

    // Original eigenvectors: v0, v1, v2

    if(eigenvectorPermutation_ == 1) {

      // Save original Q
      //Q0_ = Q;

      // Permut eigenvectors: v0, v1, v2
      //Q[0][0] = Q0_[0][0]; Q[0][1] = Q0_[0][1]; Q[0][2] = Q0_[0][2];
      //Q[1][0] = Q0_[1][0]; Q[1][1] = Q0_[1][1]; Q[1][2] = Q0_[1][2];
      //Q[2][0] = Q0_[2][0]; Q[2][1] = Q0_[2][1]; Q[2][2] = Q0_[2][2];

    } else if(eigenvectorPermutation_ == 2) {

      // Save original Q
      Q0_[0][0] = Q[0][0]; Q0_[0][1] = Q[0][1]; Q0_[0][2] = Q[0][2];
      Q0_[1][0] = Q[1][0]; Q0_[1][1] = Q[1][1]; Q0_[1][2] = Q[1][2];
      Q0_[2][0] = Q[2][0]; Q0_[2][1] = Q[2][1]; Q0_[2][2] = Q[2][2];

      // Permut eigenvectors: v0, v2, v1
      Q[0][0] = Q0_[0][0]; Q[0][1] = Q0_[0][2]; Q[0][2] = Q0_[0][1];
      Q[1][0] = Q0_[1][0]; Q[1][1] = Q0_[1][2]; Q[1][2] = Q0_[1][1];
      Q[2][0] = Q0_[2][0]; Q[2][1] = Q0_[2][2]; Q[2][2] = Q0_[2][1];

    } else if(eigenvectorPermutation_ == 3) {

      // Save original Q
      Q0_[0][0] = Q[0][0]; Q0_[0][1] = Q[0][1]; Q0_[0][2] = Q[0][2];
      Q0_[1][0] = Q[1][0]; Q0_[1][1] = Q[1][1]; Q0_[1][2] = Q[1][2];
      Q0_[2][0] = Q[2][0]; Q0_[2][1] = Q[2][1]; Q0_[2][2] = Q[2][2];

      // Permut eigenvectors: v2, v1, v0
      Q[0][0] = Q0_[0][2]; Q[0][1] = Q0_[0][1]; Q[0][2] = Q0_[0][0];
      Q[1][0] = Q0_[1][2]; Q[1][1] = Q0_[1][1]; Q[1][2] = Q0_[1][0];
      Q[2][0] = Q0_[2][2]; Q[2][1] = Q0_[2][1]; Q[2][2] = Q0_[2][0];

    }

  }

}

//--------------------------------------------------------------------------
//-------- form_perturbed_stress -------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::form_perturbed_stress(
  const double (&D)[3][3], const double (&Q)[3][3], double (&A)[3][3])
{
  // A = Q*D*QT
  double QT[3][3];
  double B[3][3];

  // compute QT
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      QT[j][i] = Q[i][j];
    }
  }
  //mat-vec, B = Q*D
  matrix_matrix_multiply(Q,D,B);

  // mat-vec, A = (Q*D)*QT = B*QT
  matrix_matrix_multiply(B,QT,A);
}

//--------------------------------------------------------------------------
//-------- matrix_matrix_multiply ------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeSolverAlgorithm::matrix_matrix_multiply(
  const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3])
{
  //C = A*B
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0;
      for (int k = 0; k < 3; ++k) {
        sum = sum + A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}

//--------------------------------------------------------------------------
//-------- van_leer ---------------------------------------------------------
//--------------------------------------------------------------------------
double
AssembleMomentumEdgeSolverAlgorithm::van_leer(
  const double &dqm,
  const double &dqp,
  const double &small)
{
  double limit = (2.0*(dqm*dqp+std::abs(dqm*dqp))) /
    ((dqm+dqp)*(dqm+dqp)+small);
  return limit;
}

} // namespace nalu
} // namespace Sierra
