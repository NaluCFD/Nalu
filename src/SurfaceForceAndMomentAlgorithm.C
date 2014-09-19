/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <SurfaceForceAndMomentAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// basic c++
#include <math.h>
#include <iostream>
#include <fstream>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// SurfaceForceAndMomentAlgorithm - post process sigma_ijnjdS and tau_wall
//                                  via lumped nodal projection (area weighed)
//                                  Driven by SurfaceForceAndMomentAlgDriver
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SurfaceForceAndMomentAlgorithm::SurfaceForceAndMomentAlgorithm(
  Realm &realm,
  stk::mesh::PartVector &partVec,
  const std::string &outputFileName,
  const int &frequency,
  const std::vector<double > &parameters,
  const bool &useShifted)
  : Algorithm(realm, partVec),
    outputFileName_(outputFileName),
    frequency_(frequency),
    parameters_(parameters),
    useShifted_(useShifted),
    includeDivU_(realm.get_divU()),
    coordinates_(NULL),
    pressure_(NULL),
    pressureForce_(NULL),
    tauWall_(NULL),
    viscosity_(NULL),
    dudx_(NULL),
    exposedAreaVec_(NULL),
    assembledArea_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  pressureForce_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "pressure_force");
  tauWall_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "tau_wall");
  // extract viscosity name
  const std::string viscName = realm_.is_turbulent()
    ? "effective_viscosity_u" : "viscosity";
  viscosity_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, viscName);
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  exposedAreaVec_ = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  assembledArea_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_area_force_moment");
  // error check on params
  const size_t nDim = meta_data.spatial_dimension();
  if ( parameters_.size() > nDim )
    throw std::runtime_error("SurfaceForce: parameter length wrong; expect nDim");

  // deal with file name and banner
  if ( Env::parallel_rank() == 0 ) {
    std::ofstream myfile;
    myfile.open(outputFileName_.c_str());
    myfile << "Time    Fx    Fy    Fz    Mx    My    Mz "<< std::endl;
    myfile.close();
  }
 }

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SurfaceForceAndMomentAlgorithm::~SurfaceForceAndMomentAlgorithm()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceForceAndMomentAlgorithm::execute()
{
  // check to see if this is a valid step to process output file
  const int timeStepCount = realm_.get_time_step_count();
  const bool processMe = (timeStepCount % frequency_) == 0 ? true : false;

  // common
  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  const int nDim = meta_data.spatial_dimension();

  // nodal fields to gather
  std::vector<double> ws_pressure;
  std::vector<double> ws_viscosity;

  // master element
  std::vector<double> ws_face_shape_function;

  const double currentTime = realm_.get_current_time();

  // local force and moment; i.e., to be assembled
  double l_force_moment[6] = {0,0,0,0,0,0};

  // work force, moment and radius; i.e., to be pushed to cross_product()
  double ws_p_force[3] = {0,0,0};
  double ws_v_force[3] = {0,0,0};
  double ws_t_force[3] = {0,0,0};
  double ws_tau[3] = {0,0,0};
  double ws_moment[3] = {0,0,0};
  double ws_radius[3] = {0,0,0};

  // will need surface normal
  double ws_normal[3] = {0,0,0};

  // centroid
  double centroid[3] = {0,0,0};
  for ( size_t k = 0; k < parameters_.size(); ++k)
    centroid[k] = parameters_[k];

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;

    // algorithm related; element
    ws_pressure.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);
    ws_face_shape_function.resize(nodesPerFace*nodesPerFace);
    
    // pointers
    double *p_pressure = &ws_pressure[0];
    double *p_viscosity = &ws_viscosity[0];
    double *p_face_shape_function = &ws_face_shape_function[0];

    // shape functions
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_face_shape_function[0]);
    else
      meFC->shape_fcn(&p_face_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // face node relations
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);

      //======================================
      // gather nodal data off of face
      //======================================
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        p_pressure[ni]    = *stk::mesh::field_data(*pressure_, node);
        p_viscosity[ni] = *stk::mesh::field_data(*viscosity_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);

      for ( int ip = 0; ip < nodesPerFace; ++ip ) {

        // offsets
        const int offSetAveraVec = ip*nDim;
        const int offSetSF_face = ip*nodesPerFace;

        // interpolate to bip
        double pBip = 0.0;
        double muBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[offSetSF_face+ic];
          pBip += r*p_pressure[ic];
          muBip += r*p_viscosity[ic];
	}

        // extract nodal fields
        stk::mesh::Entity node = face_node_rels[ip];
        const double * coord = stk::mesh::field_data(*coordinates_, node );
        const double *duidxj = stk::mesh::field_data(*dudx_, node );
	double *pressureForce = stk::mesh::field_data(*pressureForce_, node );
        double *tauWall = stk::mesh::field_data(*tauWall_, node );
        const double assembledArea = *stk::mesh::field_data(*assembledArea_, node );

        // divU and aMag
        double divU = 0.0;
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j) {
          divU += duidxj[j*nDim+j];
          aMag += areaVec[offSetAveraVec+j]*areaVec[offSetAveraVec+j];
        }
        aMag = std::sqrt(aMag);

        // normal
        for ( int i = 0; i < nDim; ++i ) {
          const double ai = areaVec[offSetAveraVec+i];
          ws_normal[i] = ai/aMag;
        }

        // load radius; assemble force -sigma_ij*njdS and compute tau_ij njDs
        for ( int i = 0; i < nDim; ++i ) {
          const double ai = areaVec[offSetAveraVec+i];
          ws_radius[i] = coord[i] - centroid[i];
          // set forces
          ws_v_force[i] = 2.0/3.0*muBip*divU*includeDivU_*ai;
          ws_p_force[i] = pBip*ai;
	  pressureForce[i] += pBip*ai;
          double dflux = 0.0;
          double tauijNj = 0.0;
          const int offSetI = nDim*i;
          for ( int j = 0; j < nDim; ++j ) {
            const int offSetTrans = nDim*j+i;
            dflux += -muBip*(duidxj[offSetI+j] + duidxj[offSetTrans])*areaVec[offSetAveraVec+j];
            tauijNj += -muBip*(duidxj[offSetI+j] + duidxj[offSetTrans])*ws_normal[j];
          }
          // accumulate viscous force and set tau for component i
          ws_v_force[i] += dflux;
          ws_tau[i] = tauijNj;
        }

        // compute total force and tangential tau
        double tauTangential = 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          ws_t_force[i] = ws_p_force[i] + ws_v_force[i];
          double tauiTangential = (1.0-ws_normal[i]*ws_normal[i])*ws_tau[i];
          for ( int j = 0; j < nDim; ++j ) {
            if ( i != j )
              tauiTangential -= ws_normal[i]*ws_normal[j]*ws_tau[j];
          }
          tauTangential += tauiTangential*tauiTangential;
        }

        // assemble nodal quantities; scaled by area for L2 lumped nodal projection
        const double fac = aMag/assembledArea;
        *tauWall += std::sqrt(tauTangential)*fac;

        cross_product(&ws_t_force[0], &ws_moment[0], &ws_radius[0]);

        // assemble force and moment
        for ( int j = 0; j < 3; ++j ) {
          l_force_moment[j] += ws_t_force[j];
          l_force_moment[j+3] += ws_moment[j];
        }

      }
    }
  }

  if ( processMe ) {
    // parallel assemble and output
    double g_force_moment[6] = {};
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    // Parallel assembly of L2
    stk::all_reduce_sum(comm, &l_force_moment[0], &g_force_moment[0], 6);

    // deal with file name and banner
    if ( Env::parallel_rank() == 0 ) {
      std::ofstream myfile;
      myfile.open(outputFileName_.c_str(), std::ios_base::app);
      myfile << currentTime << " "
        << g_force_moment[0] << " " << g_force_moment[1] << " " << g_force_moment[2] << " "
        << g_force_moment[3] << " " << g_force_moment[4] << " " << g_force_moment[5] << std::endl;
      myfile.close();
    }
  }

}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceForceAndMomentAlgorithm::pre_work()
{

  // common
  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  const int nDim = meta_data.spatial_dimension();

  //======================
  // assemble area
  //======================

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
      &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets
    = realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
          ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // face node relations
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);

      for ( int ip = 0; ip < nodesPerFace; ++ip ) {
        // offsets
        const int offSetAveraVec = ip*nDim;

        // extract nodal fields
        stk::mesh::Entity node = face_node_rels[ip];
        double *assembledArea = stk::mesh::field_data(*assembledArea_, node );

        // aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j)
          aMag += areaVec[offSetAveraVec+j]*areaVec[offSetAveraVec+j];
        aMag = std::sqrt(aMag);

        // assemble nodal quantities
        *assembledArea += aMag;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- cross_product ----------------------------------------------------
//--------------------------------------------------------------------------
void
SurfaceForceAndMomentAlgorithm::cross_product(
  double *force, double *cross, double *rad)
{
  cross[0] =   rad[1]*force[2] - rad[2]*force[1];
  cross[1] = -(rad[0]*force[2] - rad[2]*force[0]);
  cross[2] =   rad[0]*force[1] - rad[1]*force[0];
}

} // namespace nalu
} // namespace Sierra
