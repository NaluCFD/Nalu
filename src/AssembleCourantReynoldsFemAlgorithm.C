/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleCourantReynoldsFemAlgorithm.h>
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
// AssembleCourantReynoldsFemAlgorithm - Courant number calc for FEM
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleCourantReynoldsFemAlgorithm::AssembleCourantReynoldsFemAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    meshMotion_(realm_.does_mesh_move()),
    velocityRTM_(NULL),
    coordinates_(NULL),
    density_(NULL),
    viscosity_(NULL),
    elemReynolds_(NULL),
    elemCourant_(NULL)
{
  // save off data
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  const std::string viscName = (realm.is_turbulent())
     ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);

  // provide for elemental fields
  elemReynolds_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_reynolds");
  elemCourant_ = meta_data.get_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_courant");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleCourantReynoldsFemAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const double scaleFac = nDim == 3 ? 3.0 : 2.0;

  const double dt = realm_.timeIntegrator_->get_time_step();
  const double small = 1.0e-16;

  // nodal fields to gather; gather everything other than what we are assembling
  std::vector<double> ws_vrtm;
  std::vector<double> ws_coordinates;
  std::vector<double> ws_density;
  std::vector<double> ws_viscosity;

  // master element space
  std::vector<double> ws_shape_function;
  std::vector<double> ws_deriv;
  std::vector<double> ws_gUpper;
  std::vector<double> ws_gLower;

  // fixed size
  std::vector<double> ws_vrtmIp(nDim);

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // set courant/reynolds number to something small
  double maxCR[2] = {-1.0, -1.0};

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( const stk::mesh::Bucket* bucket_ptr : elem_buckets ) {
    const stk::mesh::Bucket & b = *bucket_ptr;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meFEM->nodesPerElement_;
    const int numIp = meFEM->numIntPoints_;

    // algorithm related
    ws_vrtm.resize(nodesPerElement*nDim);
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_density.resize(nodesPerElement);
    ws_viscosity.resize(nodesPerElement);
    ws_shape_function.resize(numIp*nodesPerElement);
    ws_deriv.resize(nDim*numIp*nodesPerElement);
    ws_gUpper.resize(nDim*nDim*numIp); // g^ij (covariant)
    ws_gLower.resize(nDim*nDim*numIp); // g_ij (contravariat)

    // pointers.
    double *p_vrtm = &ws_vrtm[0];
    double *p_coordinates = &ws_coordinates[0];
    double *p_density = &ws_density[0];
    double *p_viscosity = &ws_viscosity[0];

    // fixed master element evaluations
    meFEM->shape_fcn(&ws_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get elem
      stk::mesh::Entity elem = b[k];
      
      // pointers to elem data
      double * elemReynolds = stk::mesh::field_data(*elemReynolds_, b, k);
      double * elemCourant = stk::mesh::field_data(*elemCourant_, b, k);
      
      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = bulk_data.begin_nodes(elem);
      int num_nodes = bulk_data.num_nodes(elem);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        double * coords = stk::mesh::field_data(*coordinates_, node );
        double * vrtm = stk::mesh::field_data(*velocityRTM_, node );

        // gather scalars
        p_density[ni]   = *stk::mesh::field_data(densityNp1, node);
        p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j = 0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
          p_vrtm[offSet+j] = vrtm[j];
        }
      }

      // compute gij; requires a proper ws_deriv from above
      meFEM->gij(&p_coordinates[0], &ws_gUpper[0], &ws_gLower[0], &ws_deriv[0]);
      
      // compute cfl and Re at each ip; set ip max to negative
      double eReynolds = -1.0;
      double eCourant = -1.0;
      for ( int ip = 0; ip < numIp; ++ip ) {
        
        // pointer to g_ij
        const double *p_gLower = &ws_gLower[nDim*nDim*ip];

        // zero out vector
        for ( int i = 0; i < nDim; ++i ) {
          ws_vrtmIp[i] = 0.0;
        }
        
        // determine ip values of interest
        double rhoIp = 0.0;
        double muIp = 0.0;
        const int ipNpe = ip*nodesPerElement;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          // save off shape function
          const double r = ws_shape_function[ipNpe+ic];
          
          rhoIp += r*p_density[ic];
          muIp += r*p_viscosity[ic];

          for ( int j = 0; j < nDim; ++j ) {
            ws_vrtmIp[j] += r*p_vrtm[ic*nDim+j];
          }
        }
        
        double vrtmiGLowerVrtmj = 0.0;
        double glSq = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          const double vrtmi = ws_vrtmIp[i];
          for ( int j = 0; j < nDim; ++j ) {
            const double gijL = p_gLower[i*nDim+j];
            vrtmiGLowerVrtmj += vrtmi*gijL*ws_vrtmIp[j];
            glSq += gijL*gijL; 
          }
        }
        const double ipCourant = std::sqrt(vrtmiGLowerVrtmj)*dt/2.0;
        maxCR[0] = std::max(maxCR[0], ipCourant);
        
        const double diffIp = muIp/rhoIp + small;
        const double ipReynolds = 2.0*std::sqrt(vrtmiGLowerVrtmj/(diffIp*diffIp*glSq/scaleFac));
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
  
  // parallel max
  double g_maxCR[2]  = {};
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_max(comm, maxCR, g_maxCR, 2);

  // sent to realm
  realm_.maxCourant_ = g_maxCR[0];
  realm_.maxReynolds_ = g_maxCR[1];

}

} // namespace nalu
} // namespace Sierra
