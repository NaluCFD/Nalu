/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleMomentumMoengWallFunctionSolverAlgorithm.h>
#include <SolverAlgorithm.h>
#include <EquationSystem.h>
#include <LinearSystem.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>
#include <ABLProfileFunction.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <cmath>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleMomentumMoengWallFunctionSolverAlgorithm - ABL utau at wall bc
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleMomentumMoengWallFunctionSolverAlgorithm::AssembleMomentumMoengWallFunctionSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  const bool &useShifted,
  const double &gravity,
  const double &z0,
  const double &Tref)
  : SolverAlgorithm(realm, part, eqSystem),
    useShifted_(useShifted),
    z0_(z0), 
    Tref_(Tref), 
    gravity_(gravity), 
    alpha_h_(1.0), 
    beta_m_(16.0), 
    beta_h_(16.0),
    gamma_m_(5.0),
    gamma_h_(5.0),
    kappa_(realm.get_turb_model_constant(TM_kappa))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  bcVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "wall_velocity_bc");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  bcHeatFlux_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "heat_flux_bc");
  specificHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  wallNormalDistanceBip_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "wall_normal_distance_bip");
  tauWallBip_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "tau_wall_bip");
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumMoengWallFunctionSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildFaceToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumMoengWallFunctionSolverAlgorithm::execute()
{

  ABLProfileFunction *p_ABLProfFun;
  StableABLProfileFunction StableProfFun(gamma_m_, gamma_h_);
  UnstableABLProfileFunction UnstableProfFun(beta_m_, beta_h_);
  NeutralABLProfileFunction NeutralProfFun;

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS; nodesPerFace*nDim*nodesPerFace*nDim and nodesPerFace*nDim
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // bip values
  std::vector<double> uBip(nDim);
  std::vector<double> uBcBip(nDim);
  std::vector<double> unitNormal(nDim);

  // pointers to fixed values
  double *p_uBip = &uBip[0];
  double *p_uBcBip = &uBcBip[0];
  double *p_unitNormal= &unitNormal[0];

  // nodal fields to gather
  std::vector<double> ws_velocityNp1;
  std::vector<double> ws_bcVelocity;
  std::vector<double> ws_bcHeatFlux;
  std::vector<double> ws_density;
  std::vector<double> ws_specificHeat;


  // master element
  std::vector<double> ws_shape_function;
  std::vector<double> ws_face_shape_function;

  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // declare and initialize surface-averaged quantities
  std::vector<double> utan_sa;
  std::vector<double> tau_sa;
  utan_sa.resize(nDim);
  tau_sa.resize(nDim);
  for (int i=0; i<nDim; ++i) {
    utan_sa[i] = 0.0;
    tau_sa[i] = 0.0;
  }
  double utau_sa = 0.0;
  double V_sa = 0.0;
  double aTot = 0.0;
  double zp_sa = 0.0;
  double rho_sa = 0.0;
  double Cp_sa = 0.0;
  double heatFlux_sa = 0.0;

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];

    // extract master element
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal; face perspective (use with face_node_relations)
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nDim*nodesPerFace*nDim;
    const int rhsSize = nodesPerFace*nDim;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerFace);

    // algorithm related; element
    ws_velocityNp1.resize(nodesPerFace*nDim);
    ws_bcVelocity.resize(nodesPerFace*nDim);
    ws_bcHeatFlux.resize(nodesPerFace);
    ws_density.resize(nodesPerFace);
    ws_specificHeat.resize(nodesPerFace);

    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_bcVelocity = &ws_bcVelocity[0];
    double *p_bcHeatFlux = &ws_bcHeatFlux[0];
    double *p_density = &ws_density[0];
    double *p_specificHeat = &ws_specificHeat[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);

    // Loop to calculate surface-averaged velocities
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;

        // gather scalars
        p_bcHeatFlux[ni] = *stk::mesh::field_data(*bcHeatFlux_, node);
        p_density[ni]    = *stk::mesh::field_data(densityNp1, node);
        p_specificHeat[ni] = *stk::mesh::field_data(*specificHeat_, node);

        // gather velocity vectors
        double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        double * uBc = stk::mesh::field_data(*bcVelocity_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_velocityNp1[offSet+j] = uNp1[j];
          p_bcVelocity[offSet+j] = uBc[j];
        }
      }

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulk_data.begin_elements(face);
      ThrowAssert( bulk_data.num_elements(face) == 1 );

      // get element; its face ordinal number
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = bulk_data.begin_element_ordinals(face)[0];

      // pointer to face data
      const double *areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      double *wallNormalDistanceBip = stk::mesh::field_data(*wallNormalDistanceBip_, face);

      // get the relations off of element
      stk::mesh::Entity const * elem_node_rels = bulk_data.begin_nodes(element);

      // loop over face nodes
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);
        const int localFaceNode = faceIpNodeMap[ip];

        // left and right nodes; right is on the face; left is the opposing node
        stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
        stk::mesh::Entity nodeR = face_node_rels[localFaceNode];

        // extract nodal fields
        const double * coordL = stk::mesh::field_data(*coordinates_, nodeL );
        const double * coordR = stk::mesh::field_data(*coordinates_, nodeR );

        // zero out vector quantities; squeeze in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
          p_uBcBip[j] = 0.0;
          const double axj = areaVec[offSetAveraVec+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate to bip
        double heatFluxBip = 0.0;
        double rhoBip = 0.0;
        double CpBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          rhoBip += r*p_density[ic];
          CpBip += r*p_specificHeat[ic];
          heatFluxBip += r*p_bcHeatFlux[ic];
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_velocityNp1[offSetFN+j];
            p_uBcBip[j] += r*p_bcVelocity[offSetFN+j];
          }
        }

        rho_sa += rhoBip*aMag;
        heatFlux_sa += heatFlux_sa*aMag;
        Cp_sa += CpBip*aMag;

        // form unit normal
        for ( int j = 0; j < nDim; ++j ) {
          p_unitNormal[j] = areaVec[offSetAveraVec+j]/aMag;
        }

        // form unit normal and determine yp (approximated by 1/4 distance along edge)
        double ypBip = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double nj = areaVec[offSetAveraVec+j]/aMag;
          const double ej = 0.25*(coordR[j] - coordL[j]);
          ypBip += nj*ej*nj*ej;
          p_unitNormal[j] = nj;
        }
        ypBip = std::sqrt(ypBip);
        wallNormalDistanceBip[ip] = ypBip;

       // determine tangential velocity
        double uTangential = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          double uiTan = 0.0;
          double uiBcTan = 0.0;
          for ( int j = 0; j < nDim; ++j ) {
            const double ninj = p_unitNormal[i]*p_unitNormal[j];
            if ( i==j ) {
              const double om_nini = 1.0 - ninj;
              uiTan += om_nini*p_uBip[j];
              uiBcTan += om_nini*p_uBcBip[j];
            }
            else {
              uiTan -= ninj*p_uBip[j];
              uiBcTan -= ninj*p_uBcBip[j];
            }
          }
          uTangential += (uiTan-uiBcTan)*(uiTan-uiBcTan);
          utan_sa[i] += uiTan*aMag; // sum into surface-averaged tangential velocity
        }
        uTangential = std::sqrt(uTangential);
        V_sa += uTangential*aMag;
        aTot += aMag;
        zp_sa += wallNormalDistanceBip[ip]*aMag;

      } //loop over face integration points
    } // loop over faces
  } // loop over buckets 

  // parallel assemble surface-averaged quantities
  double g_utan_sa[nDim];
  for ( int i = 0; i < nDim; ++i ) {
    g_utan_sa[i] = 0.0;
  }
  double g_V_sa = 0.0;
  double g_aTot = 0.0;
  double g_zp_sa = 0.0;
  double g_rho_sa = 0.0;
  double g_heatFlux_sa = 0.0;
  double g_Cp_sa = 0.0;
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_sum(comm, &utan_sa[0], &g_utan_sa[0], nDim);
  stk::all_reduce_sum(comm, &V_sa, &g_V_sa, 1);
  stk::all_reduce_sum(comm, &aTot, &g_aTot, 1);
  stk::all_reduce_sum(comm, &zp_sa, &g_zp_sa, 1);
  stk::all_reduce_sum(comm, &rho_sa, &g_rho_sa, 1);
  stk::all_reduce_sum(comm, &heatFlux_sa, &g_heatFlux_sa, 1);
  stk::all_reduce_sum(comm, &Cp_sa, &g_Cp_sa, 1);
  for ( int i=0; i < nDim; ++i ) {
    g_utan_sa[i] /= g_aTot;
  }
  g_V_sa /= g_aTot;
  g_zp_sa /= g_aTot;
  g_rho_sa /= g_aTot;
  g_heatFlux_sa /= g_aTot;
  g_Cp_sa /= g_aTot;

  // determine surface-averaged shear stress from Monin-Obukhov expressions
  const double Tflux = g_heatFlux_sa / (g_rho_sa * g_Cp_sa);

  const double eps_heat_flux = 1.0e-8;
  const double largenum = 1.0e8;
  double Lfac;
  if (Tflux < -eps_heat_flux) {
    p_ABLProfFun = &StableProfFun;
  }
  else if (Tflux > eps_heat_flux) {
    p_ABLProfFun = &UnstableProfFun;
  }
  else {
    p_ABLProfFun = &NeutralProfFun;
  }
  if (std::abs(Tflux) < eps_heat_flux) {
    Lfac = largenum;
  }
  else {
    Lfac = -Tref_ / (kappa_ * gravity_ * Tflux);
  }

  compute_utau(g_V_sa, g_zp_sa, Tflux, p_ABLProfFun, utau_sa);
  double L = utau_sa*utau_sa*utau_sa * Lfac;
  double lambda = g_rho_sa*kappa_*utau_sa/(std::log(g_zp_sa/z0_) - p_ABLProfFun->velocity(g_zp_sa/L));
  for (int i=0; i<nDim; ++i) {
    tau_sa[i] = lambda*g_utan_sa[i];
  }

  // Loop to calculate boundary condition terms
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal; face perspective (use with face_node_relations)
    const int *faceIpNodeMap = meFC->ipNodeMap();

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nDim*nodesPerFace*nDim;
    const int rhsSize = nodesPerFace*nDim;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerFace);

    // algorithm related; element
    ws_velocityNp1.resize(nodesPerFace*nDim);

    ws_face_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_velocityNp1 = &ws_velocityNp1[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;

        // gather vectors
        double * uNp1 = stk::mesh::field_data(velocityNp1, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_velocityNp1[offSet+j] = uNp1[j];
        }
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      double *tauWallBip = stk::mesh::field_data(*tauWallBip_, face);

      // loop over face nodes
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        const int localFaceNode = faceIpNodeMap[ip];

        // zero out vector quantities; squeeze in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_uBip[j] = 0.0;
          const double axj = areaVec[offSetAveraVec+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate to bip
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          const int offSetFN = ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_uBip[j] += r*p_velocityNp1[offSetFN+j];
          }
        }

        // form unit normal
        for ( int j = 0; j < nDim; ++j ) {
          p_unitNormal[j] = areaVec[offSetAveraVec+j]/aMag;
        }

        // start the lhs assembly
        double lambda_1[nDim];
        double lambda_2[nDim];
        double uiTan[nDim];
        for ( int i = 0; i < nDim; ++i) {
          double sgnu = (g_utan_sa[i] > 0.0)?1.0:-1.0;
          lambda_1[i] = tau_sa[i] / std::max(g_V_sa,1.0e-10) * aMag;
          lambda_2[i] = tau_sa[i] / (sgnu*std::max(std::abs(g_utan_sa[i]),1.0e-10)) * aMag;
          uiTan[i] = 0.0;
        }
        double uTangential = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          for ( int j = 0; j < nDim; ++j ) {
            const double ninj = p_unitNormal[i]*p_unitNormal[j];
            if ( i==j ) {
              const double om_nini = 1.0 - ninj;
              uiTan[i] += om_nini*p_uBip[j];
            }
            else {
              uiTan[i] -= ninj*p_uBip[j];
            }
          }
          uTangential += uiTan[i]*uiTan[i];
        }   
        uTangential = std::sqrt(uTangential);

        for ( int i = 0; i < nDim; ++i ) {

          int indexR = localFaceNode*nDim + i;
          int rowR = indexR*nodesPerFace*nDim;

          lambda_1[i] /= std::max(1.0e-10,uTangential);
          for ( int j=0; j < nDim; ++j ) {
            for ( int k=0; k < nDim; ++k ) {
              const double njnk = p_unitNormal[j]*p_unitNormal[k];
              if ( j==k ) {
                const double om_njnk = 1.0 - njnk;
                for (int ic = 0; ic < nodesPerFace; ++ic)
                  p_lhs[rowR+ic*nDim+j] += lambda_1[i]*uiTan[k]*om_njnk*p_face_shape_function[offSetSF_face+ic];
              }
              else {
                for (int ic = 0; ic < nodesPerFace; ++ic)
                  p_lhs[rowR+ic*nDim+j] -= lambda_1[i]*uiTan[k]*njnk*p_face_shape_function[offSetSF_face+ic];
              }
            }
          }

          for ( int j = 0; j < nDim; ++j ) {
            const double ninj = p_unitNormal[i]*p_unitNormal[j];
            if ( i==j ) {
              const double om_nini = 1.0 - ninj;
              for (int ic = 0; ic < nodesPerFace; ++ic)
                p_lhs[rowR+ic*nDim+i] += lambda_2[i]*om_nini*p_face_shape_function[offSetSF_face+ic];
            }
            else {
              for (int ic = 0; ic < nodesPerFace; ++ic)
                p_lhs[rowR+ic*nDim+j] -= lambda_2[i]*ninj*p_face_shape_function[offSetSF_face+ic];
            }
          }
          double sgnu = (g_utan_sa[i] > 0.0)?1.0:-1.0;
          double tau_s_i = tau_sa[i] * ( (uTangential*g_utan_sa[i] + g_V_sa*(uiTan[i]-g_utan_sa[i])) / (sgnu*std::max(1.0e-10, std::abs(g_V_sa*g_utan_sa[i]))) );
          // Store tauWall for forces and moments algorithm to use
          tauWallBip[ip*nDim+i] = tau_s_i;
          p_rhs[indexR] -= tau_s_i*aMag;
        }
      }

      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_utau----------------------------------------------
//--------------------------------------------------------------------------
void
AssembleMomentumMoengWallFunctionSolverAlgorithm::compute_utau(
    const double &up, const double &zp, const double &qsurf, const ABLProfileFunction *ABLProfFun, double &utau )
{
  bool converged = false;
  const int maxIteration = 40;
  const double tolerance = 1.0e-7;

  double Lfac;
  if (qsurf == 0.0) {
    Lfac = 0.0;
  }
  else {
    Lfac = -Tref_ / (kappa_ * gravity_ * qsurf);
  }
  const double log_z_over_z0 = std::log(zp / z0_);

  // initial guesses for utau
  const double eps_u = 1.0e-8;
  const double perturb = 1.0e-3;
 
  double utau0;
  if (std::abs(up) < eps_u) {
    utau0 = eps_u;
    utau = eps_u;
    return;
  }
  else {
    utau0 = kappa_ * up / log_z_over_z0;
  }
  if (qsurf > 0.0) { // if unstable ABL
    utau0 = 3*utau0; // push initial guess above the singularity in the function to be zero'd
  }
  double utau1 = (1.0+perturb) * utau0;

  for ( int k = 0; k < maxIteration; ++k) {

    // calculate Monin-Obukhov length
    double L0 = utau0*utau0*utau0 * Lfac;
    double L1 = utau1*utau1*utau1 * Lfac;

    // limit the values of L...
    //   to be negative and finite for qsurf>0 (unstable)
    //   to be positive and finite for qsurf<0 (stable)
    double sgnq = (qsurf > 0.0)?1.0:-1.0;
    L0 = - sgnq * std::max(1.0e-10,std::abs(L0));
    L1 = - sgnq * std::max(1.0e-10,std::abs(L1));

    // calculate normalized coordinate
    const double znorm0 = zp / L0;
    const double znorm1 = zp / L1;

    const double denom0 = log_z_over_z0 - ABLProfFun->velocity(znorm0);
    const double denom1 = log_z_over_z0 - ABLProfFun->velocity(znorm1);

    // form function to be zeroed
    const double f0 = utau0 - up * kappa_ / denom0;
    const double f1 = utau1 - up * kappa_ / denom1;

    // estimate slope, d f/d utau
    double dutau = utau1 - utau0;
    if (dutau > 0.0) {
      dutau  = std::max(1.0e-15, dutau);
    }
    else  {
      dutau = std::min(-1.0e-15, dutau);
    }
    double fPrime = (f1 - f0) / dutau;
    if (fPrime > 0.0) {
      fPrime  = std::max(1.0e-15, fPrime);
    }
    else  {
      fPrime = std::min(-1.0e-15, fPrime);
    }

    // update utau
    const double utau_tmp = utau1;
    utau1 = utau0 - f0 / fPrime;
    //enforce non-negativity of utau
    //utau1 = std::max(0.0, utau1);
    utau0 = utau_tmp;

    // check for convergence
    if ( (std::abs(f1) < tolerance) ) {
      converged = true;
      //enforce non-negativity of utau
      utau = std::max(0.0, utau1);
      break;
    }
  }

  // report trouble
  if ( !converged ) {
    NaluEnv::self().naluOutputP0() << "Issue with ABL utau: not converged " << std::endl;
    NaluEnv::self().naluOutputP0() << up << " " << zp << " " << utau0 << " " << utau1 << std::endl;
  }


  // NOTE: SOWFA implementation contains another block to check if utau
  // converged to a small value.  In this case, a different zero to
  // the function is sought using a different algorithm.
}


} // namespace nalu
} // namespace Sierra
