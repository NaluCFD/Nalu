/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleCourantReynoldsElemAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <NaluEnv.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleCourantReynoldsElemAlgorithm - Courant number calc for both edge
//                                        and elem
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleCourantReynoldsElemAlgorithm::AssembleCourantReynoldsElemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part, 
  MasterElement *meSCS)
  : Algorithm(realm, part),
    meSCS_(meSCS),
    nodesPerElement_(meSCS->nodesPerElement_),
    numScsIp_(meSCS->numIntPoints_),
    lrscv_(meSCS->adjacentNodes()),
    bulkData_(&realm.bulk_data()),
    metaData_(&realm_.meta_data()),
    nDim_(realm_.meta_data().spatial_dimension()),
    meshMotion_(realm_.does_mesh_move()),
    velocityRTM_(NULL),
    coordinates_(NULL),
    density_(NULL),
    viscosity_(NULL),
    elemReynolds_(NULL),
    elemCourant_(NULL)
{
  // save off data
  if ( meshMotion_ )
    velocityRTM_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  const std::string viscName = (realm.is_turbulent())
    ? "effective_viscosity_u" : "viscosity";
  viscosity_ = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);

  // provide for elemental fields
  elemReynolds_ = metaData_->get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_reynolds");
  elemCourant_ = metaData_->get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_courant");

  // size the workset fields
  ws_vrtm_.resize(nodesPerElement_*nDim_);
  ws_coordinates_.resize(nodesPerElement_*nDim_);
  ws_density_.resize(nodesPerElement_);
  ws_viscosity_.resize(nodesPerElement_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleCourantReynoldsElemAlgorithm::execute()
{
  // pointers.
  double *p_vrtm = &ws_vrtm_[0];
  double *p_coordinates = &ws_coordinates_[0];
  double *p_density = &ws_density_[0];
  double *p_viscosity = &ws_viscosity_[0];

  const double dt = realm_.timeIntegrator_->get_time_step();
  const double small = 1.0e-16;

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // set courant/reynolds number to something small
  double maxCR[2] = {-1.0, -1.0};

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData_->locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get elem
      stk::mesh::Entity elem = b[k];
      
      // pointers to elem data
      double * elemReynolds = stk::mesh::field_data(*elemReynolds_, b, k);
      double * elemCourant = stk::mesh::field_data(*elemCourant_, b, k);
      
      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = bulkData_->begin_nodes(elem);
      int num_nodes = bulkData_->num_nodes(elem);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement_ );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        double * coords = stk::mesh::field_data(*coordinates_, node );
        double * vrtm = stk::mesh::field_data(*velocityRTM_, node );

        // gather scalars
        p_density[ni]   = *stk::mesh::field_data(densityNp1, node);
        p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);

        // gather vectors
        const int offSet = ni*nDim_;
        for ( int j = 0; j < nDim_; ++j ) {
          p_coordinates[offSet+j] = coords[j];
          p_vrtm[offSet+j] = vrtm[j];
        }
      }

      // compute cfl and Re along each edge; set ip max to negative
      double eReynolds = -1.0;
      double eCourant = -1.0;
      for ( int ip = 0; ip < numScsIp_; ++ip ) {

        // left and right nodes for this ip
        const int il = lrscv_[2*ip];
        const int ir = lrscv_[2*ip+1];

        double udotx = 0.0;
        double dxSq = 0.0;
        for ( int j = 0; j < nDim_; ++j ) {
          double ujIp = 0.5*(p_vrtm[il*nDim_+j]+p_vrtm[il*nDim_+j]);
          double dxj = p_coordinates[ir*nDim_+j] - p_coordinates[il*nDim_+j];
          udotx += dxj*ujIp;
          dxSq += dxj*dxj;
        }

        udotx = std::abs(udotx);
        const double ipCourant = std::abs(udotx*dt/dxSq);
        maxCR[0] = std::max(maxCR[0], ipCourant);

        const double diffIp = 0.5*( p_viscosity[il]/p_density[il] + p_viscosity[ir]/p_density[ir] );
        const double ipReynolds = udotx/(diffIp+small);
        maxCR[1] = std::max(maxCR[1], ipReynolds);
        
        // determine local max ip value
        eReynolds = std::max(eReynolds, ipReynolds);
        eCourant = std::max(eCourant, ipCourant);
      }
      
      // scatter
      elemReynolds[0] = eReynolds;
      elemCourant[0] = eCourant;
    }
  }

  // FIXME: Tricky now since algorithms do not operate on all of the parts...

  // parallel max
  double g_maxCR[2]  = {};
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_max(comm, maxCR, g_maxCR, 2);

  // sent to realm
  realm_.maxCourant_ = std::max(realm_.maxCourant_, g_maxCR[0]);
  realm_.maxReynolds_ = std::max(realm_.maxReynolds_, g_maxCR[1]);
}

} // namespace nalu
} // namespace Sierra
