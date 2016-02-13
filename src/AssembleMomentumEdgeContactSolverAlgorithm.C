/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleMomentumEdgeContactSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>
#include <HaloInfo.h>
#include <ContactInfo.h>
#include <ContactManager.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>

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
// AssembleMomentumEdgeContactSolverAlgorithm - add LHS/RHS for momentum
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeContactSolverAlgorithm::AssembleMomentumEdgeContactSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.has_mesh_motion()),
    includeDivU_(realm_.get_divU()),
    meshVelocity_(NULL),
    pecletFunction_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    meshVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  // extract viscosity  name
  const std::string viscName = realm_.is_turbulent()
    ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  haloMdot_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_mdot");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(velocity_->name());

  // populate fieldVec; no state
  ghostFieldVec_.push_back(dudx_);
  ghostFieldVec_.push_back(coordinates_);
  // with state
  ghostFieldVec_.push_back(&(velocity_->field_of_state(stk::mesh::StateNP1)));
  ghostFieldVec_.push_back(&(density_->field_of_state(stk::mesh::StateNP1)));
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumEdgeContactSolverAlgorithm::~AssembleMomentumEdgeContactSolverAlgorithm()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeContactSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeHaloNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumEdgeContactSolverAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

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

  // space for LHS/RHS; (nodesPerElem+1)*nDim*(nodesPerElem+1)*nDim; (nodesPerElem+1)*nDim
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // space for dui/dxj. This variable is the modifed gradient with NOC
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

  // space for interpolated right state (halo)
  double densityR;
  double viscosityR;
  std::vector<double> uNp1R(nDim);
  std::vector<double> dudxR(nDim*nDim);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
  const int sizeOfVectorField = nDim;
  const int sizeOfTensorField = nDim*nDim;
  
  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // mesh motion
  std::vector<double> vrtmL(nDim);
  std::vector<double> vrtmR(nDim);
  double * p_vrtmL = &vrtmL[0];
  double * p_vrtmR = &vrtmR[0];

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
    std::vector<double> elemNodalP(nodesPerElement);
    std::vector<double> elemNodalUnp1(nDim*nodesPerElement);
    std::vector<double> elemNodalVisc(nodesPerElement);
    std::vector<double> elemNodalRho(nodesPerElement);
    std::vector<double> elemNodalDudx(nDim*nDim*nodesPerElement);
    std::vector<double> shpfc(nodesPerElement);

    // resize some things; matrix related
    const int npePlusOne = nodesPerElement+1;
    const int lhsSize = npePlusOne*nDim*npePlusOne*nDim;
    const int rhsSize = npePlusOne*nDim;
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

        elemNodalRho[ni] = *stk::mesh::field_data(densityNp1, node);
        elemNodalVisc[ni] = *stk::mesh::field_data(*viscosity_, node);

        // load up vectors/tensor
        const double *uNp1 = stk::mesh::field_data(velocityNp1, node );
        const double *dudx = stk::mesh::field_data(*dudx_, node );
        for ( int i = 0; i < nDim; ++i ) {
          const int offSet = i*nodesPerElement + ni;
          elemNodalUnp1[offSet] = uNp1[i];

          const int rowI = i*nDim;
          const int offSetT = i*nodesPerElement*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            elemNodalDudx[offSetT+j*nodesPerElement+ni] = dudx[rowI+j];
          }
        }
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const double *coordL = stk::mesh::field_data(*coordinates_, infoObject->faceNode_);
      const double *coordR = &infoObject->haloNodalCoords_[0];

      const double *dudxL = stk::mesh::field_data(*dudx_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfTensorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalDudx[0],
        &(dudxR[0]));

      const double *uNp1L = stk::mesh::field_data(velocityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalUnp1[0],
        &(uNp1R[0]));

      const double densityL = *stk::mesh::field_data(densityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalRho[0],
        &densityR);

      const double viscosityL = *stk::mesh::field_data(*viscosity_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalVisc[0],
        &viscosityR);

      // copy to velocity relative to mesh; squeeze in extrapolated values
      for ( int i = 0; i < nDim; ++i ) {
        p_vrtmL[i] = uNp1L[i];
        p_vrtmR[i] = uNp1R[i];	
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
        const double axj = p_areaVec[j];
        const double dxj = coordR[j] - coordL[j];
        axdx += axj*dxj;
        asq += axj*axj;
        udotx += 0.5*dxj*(p_vrtmL[j] + p_vrtmR[j]);
      }

      const double inv_axdx = 1.0/axdx;

      // ip props
      const double viscIp = 0.5*(viscosityL + viscosityR);
      const double diffIp = 0.5*(viscosityL/densityL + viscosityR/densityR);

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

      // divU
      double divU = 0.0;
      for ( int j = 0; j < nDim; ++j)
        divU += p_duidxj[j*nDim+j];

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
        const double uiCds  = 0.5*(uiHatL + uiHatR);

        // total advection; pressure contribution in time term expression
        const double aflux = tmdot*(pecfac*uiUpwind + om_pecfac*uiCds);

        // diffusive flux; viscous tensor doted with area vector
        double dflux = 2.0/3.0*viscIp*divU*p_areaVec[i]*includeDivU_;
        const int offSetI = nDim*i;
        for ( int j = 0; j < nDim; ++j ) {
          const int offSetTrans = nDim*j+i;
          const double axj = p_areaVec[j];
          dflux += -viscIp*(p_duidxj[offSetI+j] + p_duidxj[offSetTrans])*axj;
        }

        // residal for total flux
        meSCS->general_shape_fcn(1, &(infoObject->isoParCoords_[0]), &shpfc[0]);

        const double tflux = aflux + dflux;
        const int indexL = i;

        // setup for LHS; row left is easy
        const int rowL = indexL * npePlusOne * nDim;
        const int rLiL = rowL+indexL;

        // total flux left
        p_rhs[indexL] -= tflux;

        // for ease of reading, scale left node by nodesPerElement
        for ( int ni = 0; ni < num_nodes; ++ni ) {

          const int indexR = i + nDim*(ni+1);
          const int rLiR = rowL+indexR;

          //==============================
          // advection first
          //==============================

          // upwind advection (includes 4th); left node
          double alhsfac = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
            + 0.5*alpha*om_pecfac*tmdot;
          p_lhs[rLiL] += alhsfac*inv_nodesPerElement;

          // upwind advection (incldues 4th); right node
          alhsfac = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
            + 0.5*alpha*om_pecfac*tmdot;
          p_lhs[rLiR] += alhsfac*shpfc[ni];

          // central; left; collect terms on alpha and alphaUpw
          alhsfac = 0.5*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);
          p_lhs[rLiL] += alhsfac*inv_nodesPerElement;
          p_lhs[rLiR] += alhsfac*shpfc[ni];
          // central; right n/a

          //==============================
          // diffusion second
          //==============================
          const double axi = p_areaVec[i];

          //diffusion; row IL
          p_lhs[rLiL] -= dlhsfac*inv_nodesPerElement;
          p_lhs[rLiR] += dlhsfac*shpfc[ni];

          // diffusion; row IR; n/a

          // more diffusion; see theory manual
          for ( int j = 0; j < nDim; ++j ) {
            const double lhsfacNS = -viscIp*axi*p_areaVec[j]*inv_axdx;

            const int colL = j;
            const int colR = j + nDim*(ni+1);

            // first left; IL,IL; IL,IR
            p_lhs[rowL + colL] -= lhsfacNS*inv_nodesPerElement;
            p_lhs[rowL + colR] += lhsfacNS*shpfc[ni];

            // now right, IR,IL; IR,IR; n/a
          }
        }
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
AssembleMomentumEdgeContactSolverAlgorithm::van_leer(
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
