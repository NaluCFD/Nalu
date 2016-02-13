/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <ContactInfo.h>
#include <ContactManager.h>
#include <FieldTypeDef.h>
#include <HaloInfo.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <TimeIntegrator.h>

#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleScalarEdgeContactSolverAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeContactSolverAlgorithm::AssembleScalarEdgeContactSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx,
  ScalarFieldType *diffFluxCoeff)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.has_mesh_motion()),
    scalarQ_(scalarQ),
    dqdx_(dqdx),
    diffFluxCoeff_(diffFluxCoeff),
    meshVelocity_(NULL),
    pecletFunction_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    meshVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  haloMdot_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_mdot");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(scalarQ_->name());

  // populate fieldVec; no state
  ghostFieldVec_.push_back(scalarQ_);
  ghostFieldVec_.push_back(dqdx_);
  ghostFieldVec_.push_back(diffFluxCoeff_);
  ghostFieldVec_.push_back(coordinates_);
  // with state
  ghostFieldVec_.push_back(&(velocity_->field_of_state(stk::mesh::StateNP1)));
  ghostFieldVec_.push_back(&(density_->field_of_state(stk::mesh::StateNP1)));
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEdgeContactSolverAlgorithm::~AssembleScalarEdgeContactSolverAlgorithm()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEdgeContactSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeHaloNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEdgeContactSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  const double small = 1.0e-16;

  // extract user advection options (allow to potentially change over time)
  const std::string dofName = scalarQ_->name();
  const double alpha = realm_.get_alpha_factor(dofName);
  const double alphaUpw = realm_.get_alpha_upw_factor(dofName);
  const double hoUpwind = realm_.get_upw_factor(dofName);
  const bool useLimiter = realm_.primitive_uses_limiter(dofName);

  // one minus flavor
  const double om_alpha = 1.0-alpha;
  const double om_alphaUpw = 1.0-alphaUpw;

  // space for LHS/RHS (nodesPerElem+1)*(nodesPerElem+1); nodesPerElem+1
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // deal with state
  ScalarFieldType &scalarQNp1  = scalarQ_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // mesh motion
  std::vector<double> vrtmL(nDim);
  std::vector<double> vrtmR(nDim);
  double * p_vrtmL = &vrtmL[0];
  double * p_vrtmR = &vrtmR[0];

  // space for interpolated right state (halo)
  double qNp1R;
  double densityR;
  double diffFluxCoeffR;
  std::vector<double> uNp1R(nDim);
  std::vector<double> dqdxR(nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;

  // parallel communicate ghosted entities
  if ( NULL != realm_.contactManager_->contactGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.contactManager_->contactGhosting_), ghostFieldVec_);
  
  // iterate contactInfoVec_
  std::vector<ContactInfo *>::iterator ii;
  for( ii=realm_.contactManager_->contactInfoVec_.begin();
       ii!=realm_.contactManager_->contactInfoVec_.end(); ++ii ) {
    
    // get master element type for this contactInfo
    MasterElement *meSCS  = (*ii)->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalQ(nodesPerElement);
    std::vector <double > elemNodalUnp1(nDim*nodesPerElement);
    std::vector <double > elemNodalRho(nodesPerElement);
    std::vector <double > elemNodalDiffFluxCoeff(nodesPerElement);
    std::vector <double > elemNodalDqdx(nDim*nodesPerElement);
    std::vector <double > shpfc(nodesPerElement);

    // resize some things; matrix related
    const int npePlusOne = nodesPerElement+1;
    const int lhsSize = npePlusOne*npePlusOne;
    const int rhsSize = npePlusOne;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(npePlusOne);

    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];

    // scaling for lhs
    const double inv_nodesPerElement = 1.0/double(nodesPerElement);

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // zeroing of lhs/rhs
      for ( int k = 0; k < lhsSize; ++k ) {
        p_lhs[k] = 0.0;
      }
      for ( int k = 0; k < rhsSize; ++k ) {
        p_rhs[k] = 0.0;
      }

      // pointer to edge area vector and mdot
      const double *p_areaVec = &infoObject->haloEdgeAreaVec_[0];
      const double tmdot = *stk::mesh::field_data(*haloMdot_, infoObject->faceNode_);

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem  = infoObject->owningElement_;

      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(elem);
      const int num_nodes = bulk_data.num_nodes(elem);

      // now load the elemental values for future interpolation; fill in connected nodes
      connected_nodes[0] = infoObject->faceNode_;
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        connected_nodes[ni+1] = node;

        elemNodalQ[ni] = *stk::mesh::field_data(scalarQNp1, node);
        elemNodalRho[ni] = *stk::mesh::field_data(densityNp1, node);
        elemNodalDiffFluxCoeff[ni] = *stk::mesh::field_data(*diffFluxCoeff_, node);

        // load up vectors
        const double * uNp1 = stk::mesh::field_data(velocityNp1, node );
        const double * Gjq = stk::mesh::field_data(*dqdx_, node );
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          elemNodalDqdx[offSet] = Gjq[j];
          elemNodalUnp1[offSet] = uNp1[j];
        }
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const double *coordL = stk::mesh::field_data(*coordinates_, infoObject->faceNode_);
      const double *coordR = &infoObject->haloNodalCoords_[0];

      const double *dqdxL = stk::mesh::field_data(*dqdx_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalDqdx[0],
        &(dqdxR[0]));

      const double *uNp1L = stk::mesh::field_data(velocityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalUnp1[0],
        &(uNp1R[0]));

      const double qNp1L = *stk::mesh::field_data( scalarQNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalQ[0],
        &qNp1R);

      const double densityL = *stk::mesh::field_data(densityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalRho[0],
        &densityR);

      const double diffFluxCoeffL = *stk::mesh::field_data(*diffFluxCoeff_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalDiffFluxCoeff[0],
        &diffFluxCoeffR);

      // copy to velocity relative to mesh
      for ( int j = 0; j < nDim; ++j ) {
        p_vrtmL[j] = uNp1L[j];
        p_vrtmR[j] = uNp1R[j];
      }

      // deal with mesh motion
      if ( meshMotion_ ) {
        const double * meshVelocityL = stk::mesh::field_data(*meshVelocity_, infoObject->faceNode_);
        const double * meshVelocityR = &(infoObject->haloMeshVelocity_[0]);
        for (int j = 0; j < nDim; ++j ) {
          p_vrtmL[j] -= meshVelocityL[j];
          p_vrtmR[j] -= meshVelocityR[j];
        }
      }

      // compute geometry
      double axdx = 0.0;
      double asq = 0.0;
      double udotx = 0.0;
      for (int j = 0; j < nDim; ++j ) {
        const double dxj = coordR[j] - coordL[j];
        const double axj = p_areaVec[j];
        axdx += axj*dxj;
        asq += axj*axj;
	udotx += 0.5*dxj*(p_vrtmL[j] + p_vrtmR[j]);
      }

      const double inv_axdx = 1.0/axdx;

      // ip props
      const double viscIp = 0.5*(diffFluxCoeffL + diffFluxCoeffR);
      const double diffIp = 0.5*(diffFluxCoeffL/densityL + diffFluxCoeffR/densityR);

      // Peclet factor
      const double pecfac = pecletFunction_->execute(std::abs(udotx)/(diffIp+small));
      const double om_pecfac = 1.0-pecfac;

      // left and right extrapolation; add in diffusion calc
      double dqL = 0.0;
      double dqR = 0.0;
      double nonOrth = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double dxj = coordR[j] - coordL[j];
        dqL += 0.5*dxj*dqdxL[j];
        dqR += 0.5*dxj*dqdxR[j];

        // now non-orth (over-relaxed procedure of Jasek)
        const double axj = p_areaVec[j];
        const double kxj = axj - asq*inv_axdx*dxj;
        const double GjIp = 0.5*(dqdxL[j] + dqdxR[j]);
        nonOrth += -viscIp*kxj*GjIp;
      }

      // add limiter if appropriate
      double limitL = 1.0;
      double limitR = 1.0;
      if ( useLimiter ) {
        const double dq = qNp1R - qNp1L;
        const double dqMl = 2.0*2.0*dqL - dq;
        const double dqMr = 2.0*2.0*dqR - dq;
        limitL = van_leer(dqMl, dq, small);
        limitR = van_leer(dqMr, dq, small);
      }

      // extrapolated; for now limit
      const double qIpL = qNp1L + dqL*hoUpwind*limitL;
      const double qIpR = qNp1R - dqR*hoUpwind*limitR;

      meSCS->general_shape_fcn(1, &(infoObject->isoParCoords_[0]), &shpfc[0]);

      //====================================
      // diffusive flux
      //====================================
      double lhsfac = -viscIp*asq*inv_axdx;
      double diffFlux = lhsfac*(qNp1R - qNp1L) + nonOrth;

      // left node (face)
      p_lhs[0] -=  lhsfac;
      p_rhs[0] -=  diffFlux;

      // apply IL,IR
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        p_lhs[ni+1] += lhsfac*shpfc[ni];
      }

      //====================================
      // advective flux
      //====================================

      // 2nd order central
      const double qIp = 0.5*( qNp1L + qNp1R );

      // upwind
      const double qUpwind = (tmdot > 0) ? alphaUpw*qIpL + om_alphaUpw*qIp
          : alphaUpw*qIpR + om_alphaUpw*qIp;

      // generalized central (2nd and 4th order)
      const double qHatL = alpha*qIpL + om_alpha*qIp;
      const double qHatR = alpha*qIpR + om_alpha*qIp;
      const double qCds = 0.5*(qHatL + qHatR);

      // total advection; pressure contribution in time term expression
      const double aflux = tmdot*(pecfac*qUpwind + om_pecfac*qCds);

      // left node (face)
      p_rhs[0] -= aflux;

      // for ease of reading, scale left node by nodesPerElement
      for ( int ni = 0; ni < num_nodes; ++ni ) {

        // upwind advection (includes 4th); left node
        double alhsfac = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[0] += alhsfac*inv_nodesPerElement;

        // upwind advection; right node
        alhsfac = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
          + 0.5*alpha*om_pecfac*tmdot;
        p_lhs[ni+1] += alhsfac*shpfc[ni];

        // central; left; collect terms on alpha; alphaUpw and beta
        alhsfac = 0.5*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);
        p_lhs[0] += alhsfac*inv_nodesPerElement;
        p_lhs[ni+1] += alhsfac*shpfc[ni];
        // central; right n/a
      }

      // apply to linear system
      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
    }
  }
}

//--------------------------------------------------------------------------
//-------- van_leer ---------------------------------------------------------
//--------------------------------------------------------------------------
double
AssembleScalarEdgeContactSolverAlgorithm::van_leer(
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
