/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeMdotVofElemAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ComputeMdotVofElemAlgorithm - interior mdot for bf elem continuity
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMdotVofElemAlgorithm::ComputeMdotVofElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  const SolutionOptions &solnOpts)
  : Algorithm(realm, part),
    meshMotion_(solnOpts.does_mesh_move()),
    velocityRTM_(NULL),
    Gpdx_(NULL),
    coordinates_(NULL),
    pressure_(NULL),
    density_(NULL),
    interfaceCurvature_(NULL),
    surfaceTension_(NULL),
    vof_(NULL),
    massFlowRate_(NULL),
    volumeFlowRate_(NULL),
    shiftMdot_(solnOpts.cvfemShiftMdot_),
    shiftPoisson_(solnOpts.get_shifted_grad_op("pressure")),
    buoyancyWeight_(solnOpts.buoyancyPressureStab_ ? 1.0 : 0.0)
{
   // extract fields; nodal
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "velocity");
  Gpdx_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dpdx");
  coordinates_ = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  pressure_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "pressure");
  density_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  interfaceCurvature_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "interface_curvature");
  surfaceTension_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "surface_tension");
  vof_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  massFlowRate_ = meta_data.get_field<double>(stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");
  volumeFlowRate_ = meta_data.get_field<double>(stk::topology::ELEMENT_RANK, "volume_flow_rate_scs");
  gravity_ = realm_.solutionOptions_->gravity_;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ComputeMdotVofElemAlgorithm::~ComputeMdotVofElemAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMdotVofElemAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;
  
  // nodal fields to gather
  std::vector<double> ws_vrtm;
  std::vector<double> ws_Gpdx;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_pressure;
  std::vector<double> ws_density;
  std::vector<double> ws_kappa;
  std::vector<double> ws_sigma;
  std::vector<double> ws_vof;

  // geometry related to populate
  std::vector<double> ws_scs_areav;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;
  std::vector<double> ws_shape_function;

  // integration point data that depends on size
  std::vector<double> uIp(nDim);
  std::vector<double> GpdxIp(nDim);
  std::vector<double> dpdxIp(nDim);

  // pointers to everyone...
  double *p_uIp = &uIp[0];
  double *p_GpdxIp = &GpdxIp[0];
  double *p_dpdxIp = &dpdxIp[0];

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_)  
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;

    // algorithm related
    ws_vrtm.resize(nodesPerElement*nDim);
    ws_Gpdx.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_pressure.resize(nodesPerElement);
    ws_density.resize(nodesPerElement);
    ws_kappa.resize(nodesPerElement);
    ws_sigma.resize(nodesPerElement);
    ws_vof.resize(nodesPerElement);
    ws_scs_areav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_det_j.resize(numScsIp);
    ws_shape_function.resize(numScsIp*nodesPerElement);

    // pointers
    double *p_vrtm = &ws_vrtm[0];
    double *p_Gpdx = &ws_Gpdx[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_pressure = &ws_pressure[0];
    double *p_density = &ws_density[0];
    double *p_kappa = &ws_kappa[0];
    double *p_sigma = &ws_sigma[0];
    double *p_vof = &ws_vof[0];
    double *p_scs_areav = &ws_scs_areav[0];
    double *p_dndx = &ws_dndx[0];
    double *p_shape_function = &ws_shape_function[0];
    
    if ( shiftMdot_)
      meSCS->shifted_shape_fcn(&p_shape_function[0]);
    else
      meSCS->shape_fcn(&p_shape_function[0]);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // pointers to elem data
      double * mdot = stk::mesh::field_data(*massFlowRate_, b, k );
      double * vdot = stk::mesh::field_data(*volumeFlowRate_, b, k );

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      STK_ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * vrtm   = stk::mesh::field_data(*velocityRTM_, node);
        const double * Gjp    = stk::mesh::field_data(*Gpdx_, node);
        const double * coords = stk::mesh::field_data(*coordinates_, node);

        // gather scalars
        p_pressure[ni] = *stk::mesh::field_data(*pressure_, node);
        p_density[ni]  = *stk::mesh::field_data(densityNp1, node);
        p_kappa[ni]  = *stk::mesh::field_data(*interfaceCurvature_, node);
        p_sigma[ni]  = *stk::mesh::field_data(*surfaceTension_, node);
        p_vof[ni]  = *stk::mesh::field_data(*vof_, node);

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_vrtm[offSet+j] = vrtm[j];
          p_Gpdx[offSet+j] = Gjp[j];
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // compute geometry
      double scs_error = 0.0;
      meSCS->determinant(1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

      // compute dndx
      if (shiftPoisson_)
        meSCS->shifted_grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);
      else
        meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scs_error);
      
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        // setup for ip values
        double rhoIp = 0.0;
        double dvofdaIp = 0.0;
        double sigmaKappaIp = 0.0;

        for ( int j = 0; j < nDim; ++j ) {
          p_uIp[j] = 0.0;
          p_GpdxIp[j] = 0.0;
          p_dpdxIp[j] = 0.0;
        }

        const int offSet = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {

          const double r = p_shape_function[offSet+ic];
          rhoIp += r*p_density[ic];
          sigmaKappaIp += r*p_sigma[ic]*p_kappa[ic];

          const double pressureIc = p_pressure[ic];
          const double vofIc = p_vof[ic];

          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            p_GpdxIp[j] += r*p_Gpdx[nDim*ic+j];
            p_uIp[j] += r*p_vrtm[nDim*ic+j];
            p_dpdxIp[j] += p_dndx[offSetDnDx+j]*pressureIc;
            dvofdaIp += p_dndx[offSetDnDx+j]*vofIc*p_scs_areav[ip*nDim+j];
          }
        }

        // assemble flow rate
        double tvdot = projTimeScale*sigmaKappaIp*dvofdaIp/rhoIp;
        for ( int j = 0; j < nDim; ++j ) {
          // balanced force approach
          tvdot += (p_uIp[j] - projTimeScale*((p_dpdxIp[j]-buoyancyWeight_*rhoIp*gravity_[j])/rhoIp - p_GpdxIp[j]))*p_scs_areav[ip*nDim+j];
        }

        // save off scs ip flow rates
        mdot[ip] = rhoIp*tvdot;
        vdot[ip] = tvdot;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
