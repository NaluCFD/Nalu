/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleScalarEigenEdgeSolverAlgorithm.h>
#include <Enums.h>
#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SolutionOptions.h>

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
// AssembleScalarEigenEdgeSolverAlgorithm - add LHS/RHS for scalar using 
//                                          GGDH with eigenvalue perturbation
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEigenEdgeSolverAlgorithm::AssembleScalarEigenEdgeSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx,
  ScalarFieldType *thermalCond,
  ScalarFieldType *specHeat,
  ScalarFieldType *turbViscosity,
  const double turbSigma)
  : SolverAlgorithm(realm, part, eqSystem),
    meshMotion_(realm_.does_mesh_move()),
    includeDivU_(realm_.get_divU()),
    turbSigma_(turbSigma),
    scalarQ_(scalarQ),
    dqdx_(dqdx),
    thermalCond_(thermalCond),
    specHeat_(specHeat),
    turbViscosity_(turbViscosity),
    velocityRTM_(NULL),
    coordinates_(NULL),
    density_(NULL),
    massFlowRate_(NULL),
    edgeAreaVec_(NULL),
    turbKe_(NULL),
    velocity_(NULL),
    dudx_(NULL),
    pecletFunction_(NULL),
    cGGDH_(3.0/2.0*realm_.get_turb_model_constant(TM_cMu)/turbSigma),
    deltaB_(realm_.solutionOptions_->eigenvaluePerturbDelta_),
    perturbTurbKe_(realm_.solutionOptions_->eigenvaluePerturbTurbKe_)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  massFlowRate_ = meta_data.get_field<ScalarFieldType>(stk::topology::EDGE_RANK, "mass_flow_rate");
  edgeAreaVec_ = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");

  // EXTRA GGDH
  turbKe_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");

  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function<double>(scalarQ_->name());

  // initialize xic
  const int biasTowards = realm_.solutionOptions_->eigenvaluePerturbBiasTowards_;
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
  else {
    BinvXt_[0] = 0.0;
    BinvXt_[1] = 0.0;
    BinvXt_[2] = 0.0;
  }

  NaluEnv::self().naluOutputP0() << "Perturbation model active: towards/delta/tke: " 
                                 << biasTowards << "/" << deltaB_ 
                                 << "/" << perturbTurbKe_ << std::endl;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarEigenEdgeSolverAlgorithm::~AssembleScalarEigenEdgeSolverAlgorithm()
{
  delete pecletFunction_;
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEigenEdgeSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildEdgeToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEigenEdgeSolverAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

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

  // space for LHS/RHS; always edge connectivity
  const int nodesPerEdge = 2;
  const int lhsSize = nodesPerEdge*nodesPerEdge;
  const int rhsSize = nodesPerEdge;
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

  // deal with state
  ScalarFieldType &scalarQNp1  = scalarQ_->field_of_state(stk::mesh::StateNP1);
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

      // get edge
      stk::mesh::Entity edge = b[k];

      stk::mesh::Entity const * edge_node_rels = bulk_data.begin_nodes(edge);

      // sanity check on number or nodes
      ThrowAssert( bulk_data.num_nodes(edge) == 2 );

      // pointer to edge area vector
      for ( int j = 0; j < nDim; ++j )
        p_areaVec[j] = av[k*nDim+j];
      const double tmdot = mdot[k];

      // left and right nodes
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];

      connected_nodes[0] = nodeL;
      connected_nodes[1] = nodeR;

      // extract nodal fields
      const double * coordL = stk::mesh::field_data(*coordinates_, nodeL);
      const double * coordR = stk::mesh::field_data(*coordinates_, nodeR);

      const double * dqdxL = stk::mesh::field_data(*dqdx_, nodeL);
      const double * dqdxR = stk::mesh::field_data(*dqdx_, nodeR);

      const double * vrtmL = stk::mesh::field_data(*velocityRTM_, nodeL);
      const double * vrtmR = stk::mesh::field_data(*velocityRTM_, nodeR);

      const double qNp1L = *stk::mesh::field_data(scalarQNp1, nodeL);
      const double qNp1R = *stk::mesh::field_data(scalarQNp1, nodeR);

      const double densityL = *stk::mesh::field_data(densityNp1, nodeL);
      const double densityR = *stk::mesh::field_data(densityNp1, nodeR);

      const double thermalCondL = *stk::mesh::field_data(*thermalCond_, nodeL);
      const double thermalCondR = *stk::mesh::field_data(*thermalCond_, nodeR);

      const double specHeatL = *stk::mesh::field_data(*specHeat_, nodeL);
      const double specHeatR = *stk::mesh::field_data(*specHeat_, nodeR);

      // EXTRA GGDH: extract projected nodal velocity gradients and velocity
      const double * dudxL = stk::mesh::field_data(*dudx_, nodeL);
      const double * dudxR = stk::mesh::field_data(*dudx_, nodeR);

      const double * uNp1L = stk::mesh::field_data(*velocity_, nodeL);
      const double * uNp1R = stk::mesh::field_data(*velocity_, nodeR);

      const double turbKeL = std::max(*stk::mesh::field_data(*turbKe_, nodeL), 1.0e-16);
      const double turbKeR = std::max(*stk::mesh::field_data(*turbKe_, nodeR), 1.0e-16);

      const double turbViscL = *stk::mesh::field_data(*turbViscosity_, nodeL);
      const double turbViscR = *stk::mesh::field_data(*turbViscosity_, nodeR);
      
      // compute geometry
      double axdx = 0.0;
      double asq = 0.0;
      double udotx = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = p_areaVec[j];
        const double dxj = coordR[j] - coordL[j];
        asq += axj*axj;
        axdx += axj*dxj;
        udotx += 0.5*dxj*(vrtmL[j] + vrtmR[j]);
      }
      const double inv_axdx = 1.0/axdx;

      // ip props
      const double rhoIp = 0.5*(densityL + densityR);
      const double lamEffectiveViscIp = 0.5*(thermalCondL/specHeatL + thermalCondR/specHeatR);
      const double nuIp = lamEffectiveViscIp/rhoIp;
      const double turbViscIp = 0.5*(turbViscL + turbViscR)/turbSigma_;
      const double turbNuIp = turbViscIp/rhoIp;
      const double turbKeIp = 0.5*(turbKeL + turbKeR);

      // EXTRA GGDH: compute duidxj
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
          duidxj_[i][j] = GjUi + (uidiff - GlUidxl)*axj*inv_axdx;
        }
      }

      // divU
      double divU = 0.0;
      for ( int j = 0; j < nDim; ++j)
        divU += duidxj_[j][j];

      // estimate a time scale
      double sijMag = 0.0;
      for ( int i = 0; i < nDim; ++i ) {
        for ( int j = 0; j < nDim; ++j ) {
          const double rateOfStrain = 0.5*(duidxj_[i][j] + duidxj_[j][i]);
          sijMag += rateOfStrain*rateOfStrain;
        }
      }
      sijMag = std::sqrt(2.0*sijMag);
      const double timeScaleIp = 1.0/sijMag;

      // now compute dqdxj
      const double qDiff = qNp1R - qNp1L;

      // start sum for NOC contribution
      double Glqdxl = 0.0;
      for ( int l = 0; l< nDim; ++l ) {
        const double dxl = coordR[l] - coordL[l];
        const double Glq = 0.5*(dqdxL[l] + dqdxR[l]);
        Glqdxl += Glq*dxl;
      }

      // form scalar gradients with NOC
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = p_areaVec[j];
        const double Gjq = 0.5*(dqdxL[j] + dqdxR[j]);
        dqdxj_[j] = Gjq + (qDiff - Glqdxl)*axj*inv_axdx;
      }

      // compute the normalized Reynolds stress
      for ( int i = 0; i < nDim; ++i ) {
        for ( int j = 0; j < nDim; ++j ) {
          const double divUTerm = ( i == j ) ? 2.0/3.0*divU*includeDivU_ : 0.0;
          b_[i][j] = (-turbNuIp*(duidxj_[i][j] + duidxj_[j][i] - divUTerm))/(2.0*turbKeIp);
        }
      }
      
      // perform the decomposition
      diagonalize(b_, Q_, D_);

      // sort D
      sort(D_);

      // perturb
      perturb(D_);

      // form new stress
      form_perturbed_stress(D_, Q_, b_);

      // remove normalization; add in tke (possibly perturbed)
      const double turbKeIpPert = std::max(turbKeIp*(1.0 + perturbTurbKe_), 1.0e-16);
      for ( int i = 0; i < nDim; ++i ) {
        for ( int j = 0; j < nDim; ++j ) {
          const double fac = ( i == j ) ? 1.0/3.0 : 0.0;
          R_[i][j] = (b_[i][j] + fac)*2.0*turbKeIpPert;
        }
      }
      // R is now the perturbed Reynolds stress; rho*R is ready to be used in flux

      // Peclet factor
      const double pecfac = pecletFunction_->execute(std::abs(udotx)/(nuIp+small));
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
        nonOrth += -lamEffectiveViscIp*kxj*GjIp;
      }

      // add limiter if appropriate
      double limitL = 1.0;
      double limitR = 1.0;
      const double dq = qNp1R - qNp1L;
      if ( useLimiter ) {
        const double dqMl = 2.0*2.0*dqL - dq;
        const double dqMr = 2.0*2.0*dqR - dq;
        limitL = van_leer(dqMl, dq, small);
        limitR = van_leer(dqMr, dq, small);
      }
      
      // extrapolated; for now limit
      const double qIpL = qNp1L + dqL*hoUpwind*limitL;
      const double qIpR = qNp1R - dqR*hoUpwind*limitR;

      //====================================
      // diffusive flux; lam lhs/rhs 
      //====================================
      double lhsfac = -lamEffectiveViscIp*asq*inv_axdx;
      double diffFlux = lhsfac*(qNp1R - qNp1L) + nonOrth;

      // add in GGDH lhs/rhs
      for ( int i = 0; i < nDim; ++i ) {
        for ( int j = 0; j < nDim; ++j ) {
          const double ggFac = -cGGDH_*timeScaleIp*rhoIp*R_[i][j]*p_areaVec[j];
          diffFlux += ggFac*dqdxj_[i];
          lhsfac += ggFac*p_areaVec[i]*inv_axdx;
        }
      }
      
      // first left
      p_lhs[0] = -lhsfac;
      p_lhs[1] = +lhsfac;
      p_rhs[0] = -diffFlux;

      // now right
      p_lhs[2] = +lhsfac;
      p_lhs[3] = -lhsfac;
      p_rhs[1] = diffFlux;

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

      // total advection
      const double aflux = tmdot*(pecfac*qUpwind + om_pecfac*qCds);

      // upwind advection (includes 4th); left node
      double alhsfac = 0.5*(tmdot+std::abs(tmdot))*pecfac*alphaUpw
        + 0.5*alpha*om_pecfac*tmdot;
      p_lhs[0] += alhsfac;
      p_lhs[2] -= alhsfac;

      // upwind advection; right node
      alhsfac = 0.5*(tmdot-std::abs(tmdot))*pecfac*alphaUpw
        + 0.5*alpha*om_pecfac*tmdot;
      p_lhs[3] -= alhsfac;
      p_lhs[1] += alhsfac;

      // central; left; collect terms on alpha and alphaUpw
      alhsfac = 0.5*tmdot*(pecfac*om_alphaUpw + om_pecfac*om_alpha);
      p_lhs[0] += alhsfac;
      p_lhs[1] += alhsfac;
      // central; right; collect terms on alpha and alphaUpw
      p_lhs[2] -= alhsfac;
      p_lhs[3] -= alhsfac;

      // total flux left
      p_rhs[0] -= aflux;
      // total flux right
      p_rhs[1] += aflux;

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }
  }
}

//--------------------------------------------------------------------------
//-------- van_leer --------------------------------------------------------
//--------------------------------------------------------------------------
double
AssembleScalarEigenEdgeSolverAlgorithm::van_leer(
  const double &dqm,
  const double &dqp,
  const double &small)
{
  double limit = (2.0*(dqm*dqp+std::abs(dqm*dqp))) /
    ((dqm+dqp)*(dqm+dqp)+small);
  return limit;
}

//--------------------------------------------------------------------------
//-------- diagonalize -----------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEigenEdgeSolverAlgorithm::diagonalize(
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
AssembleScalarEigenEdgeSolverAlgorithm::sort(
  const double (&D)[3][3])
{
  // Goal: sort diagonalization eigenvalues from high to low; save off row in D from high to low

  // fist, fill in data to sort (using a trivial N^2 "selection sort; N ~small)
  double data[3] = {D[0][0], D[1][1], D[2][2]};
  rowMap_[0] = 0;
  rowMap_[1] = 1;
  rowMap_[2] = 2;

  int j = 0;
  double tmp = 0;
  for(int i=0; i < 3; ++i){
    j = i;
    for(int k = i; k < 3; ++k){
      if(data[j] < data[k]){
        j = k;
      }
    }
    tmp = data[i];
    int tmpI = rowMap_[i];
    // deal with data
    data[i] = data[j];
    data[j] = tmp;
    // now row mapping
    rowMap_[i] = rowMap_[j];
    rowMap_[j] = tmpI;
  }
}

//--------------------------------------------------------------------------
//-------- perturb ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEigenEdgeSolverAlgorithm::perturb(
  double (&D)[3][3])
{
  // extract sorted by size (L1 > L2 > L3 )
  const double Lamdba1 = D[rowMap_[0]][rowMap_[0]];
  const double Lamdba2 = D[rowMap_[1]][rowMap_[1]];
  const double Lamdba3 = D[rowMap_[2]][rowMap_[2]];

  const double pLamdba1 = (1.0-deltaB_)*Lamdba1 + deltaB_*BinvXt_[0];
  const double pLamdba2 = (1.0-deltaB_)*Lamdba2 + deltaB_*BinvXt_[1];
  const double pLamdba3 = (1.0-deltaB_)*Lamdba3 + deltaB_*BinvXt_[2];

  D[rowMap_[0]][rowMap_[0]] = pLamdba1;
  D[rowMap_[1]][rowMap_[1]] = pLamdba2;
  D[rowMap_[2]][rowMap_[2]] = pLamdba3;
}

//--------------------------------------------------------------------------
//-------- form_perturbed_stress -------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarEigenEdgeSolverAlgorithm::form_perturbed_stress(
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
AssembleScalarEigenEdgeSolverAlgorithm::matrix_matrix_multiply(
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

} // namespace nalu
} // namespace Sierra
