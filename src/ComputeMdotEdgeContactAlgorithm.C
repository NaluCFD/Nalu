/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeMdotEdgeContactAlgorithm.h>
#include <SolverAlgorithm.h>

#include <ContactInfo.h>
#include <ContactManager.h>
#include <FieldTypeDef.h>
#include <HaloInfo.h>
#include <LinearSystem.h>
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
// ComputeMdotEdgeContactAlgorithm - compute mdot for halo alg
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeMdotEdgeContactAlgorithm::ComputeMdotEdgeContactAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    meshMotion_(realm_.has_mesh_motion()),
    meshVelocity_(NULL),
    velocity_(NULL),
    Gpdx_(NULL),
    coordinates_(NULL),
    pressure_(NULL),
    density_(NULL),
    haloMdot_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    meshVelocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  Gpdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  haloMdot_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_mdot");

  // populate fieldVec; no state
  ghostFieldVec_.push_back(Gpdx_);
  ghostFieldVec_.push_back(coordinates_);
  ghostFieldVec_.push_back(pressure_);
  // with state
  ghostFieldVec_.push_back(&(velocity_->field_of_state(stk::mesh::StateNP1)));
  ghostFieldVec_.push_back(&(density_->field_of_state(stk::mesh::StateNP1)));

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeMdotEdgeContactAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // extract noc
  const std::string dofName = "pressure";
  const double nocFac
    = (realm_.get_noc_usage(dofName) == true) ? 1.0 : 0.0;

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;

  // deal with interpolation procedure
  const double interpTogether = realm_.get_mdot_interp();
  const double om_interpTogether = 1.0-interpTogether;
  
  // mesh motion
  std::vector<double> vrtmL(nDim);
  std::vector<double> vrtmR(nDim);
  double * p_vrtmL = &vrtmL[0];
  double * p_vrtmR = &vrtmR[0];  
  
  // deal with state
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // space for interpolated right state (halo)
  double pressureR;
  double densityR;
  std::vector<double> GpdxR(nDim);
  std::vector<double> velocityNp1R(nDim);

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
    std::vector <double > elemNodalP(nodesPerElement);
    std::vector <double > elemNodalUnp1(nDim*nodesPerElement);
    std::vector <double > elemNodalRho(nodesPerElement);
    std::vector <double > elemNodalGpdx(nDim*nodesPerElement);
    std::vector <double > shpfc(nodesPerElement);

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem  = infoObject->owningElement_;

      // pointer to edge area vector
      const double *p_areaVec = &infoObject->haloEdgeAreaVec_[0];

      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(elem);
      const int num_nodes = bulk_data.num_nodes(elem);

      // now load the elemental values for future interpolation
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        elemNodalP[ni] = *stk::mesh::field_data(*pressure_, node);
        elemNodalRho[ni] = *stk::mesh::field_data(densityNp1, node);

        // load up vectors
        const double * uNp1 = stk::mesh::field_data(velocityNp1, node );
        const double * Gjp = stk::mesh::field_data(*Gpdx_, node );
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          elemNodalGpdx[offSet] = Gjp[j];
          elemNodalUnp1[offSet] = uNp1[j];
        }
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const double *coordL = stk::mesh::field_data(*coordinates_, infoObject->faceNode_);
      const double *coordR = &infoObject->haloNodalCoords_[0];

      // interpolate nodal values to point-in-elem
      const double pressureL = *stk::mesh::field_data(*pressure_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalP[0],
        &pressureR);

      const double densityL = *stk::mesh::field_data(densityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalRho[0],
        &densityR);

      const double *velocityNp1L = stk::mesh::field_data(velocityNp1, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalUnp1[0],
        &(velocityNp1R[0]));

      const double *GpdxL = stk::mesh::field_data(*Gpdx_, infoObject->faceNode_);
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalGpdx[0],
        &(GpdxR[0]));
      
      // copy to velocity relative to mesh
      for ( int j = 0; j < nDim; ++j ) {
        p_vrtmL[j] = velocityNp1L[j];
        p_vrtmR[j] = velocityNp1R[j];
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
      for (int j = 0; j < nDim; ++j ) {
        const double dxj = coordR[j] - coordL[j];
        const double axj = p_areaVec[j];
        axdx += axj*dxj;
        asq += axj*axj;
      }

      const double inv_axdx = 1.0/axdx;
      const double rhoIp = 0.5*(densityR + densityL);
      
      // mdot
      double tmdot = -projTimeScale*(pressureR - pressureL)*asq*inv_axdx;
      for ( int j = 0; j < nDim; ++j ) {
        const double axj = p_areaVec[j];
        const double dxj = coordR[j] - coordL[j];
        const double kxj = axj - asq*inv_axdx*dxj; // NOC
        const double rhoUjIp = 0.5*(densityR*p_vrtmR[j] + densityL*p_vrtmL[j]);
        const double ujIp = 0.5*(p_vrtmR[j] + p_vrtmL[j]);
        const double GjIp = 0.5*(GpdxR[j] + GpdxL[j]);
        tmdot += (interpTogether*rhoUjIp + om_interpTogether*rhoIp*ujIp + projTimeScale*GjIp)*axj 
          - projTimeScale*kxj*GjIp*nocFac;
      }

      // scatter to mdot
      double *haloMdot = stk::mesh::field_data(*haloMdot_, infoObject->faceNode_);
      *haloMdot = tmdot;
    }
  }
}

} // namespace nalu
} // namespace Sierra
