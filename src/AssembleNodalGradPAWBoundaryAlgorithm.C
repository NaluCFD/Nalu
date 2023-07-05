/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradPAWBoundaryAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNodalGradPAWBoundaryAlgorithm - balanced force Gjph boundary
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradPAWBoundaryAlgorithm::AssembleNodalGradPAWBoundaryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *pressure,
  VectorFieldType *dpdx,
  const std::string bcPressureName,
  const bool overrideFacePressure)
  : Algorithm(realm, part),
    pressure_(pressure),
    dpdx_(dpdx),
    coordinates_(nullptr),
    density_(nullptr),    
    interfaceCurvature_(NULL),
    surfaceTension_(NULL),
    vof_(NULL),
    bcPressure_(nullptr),
    areaWeight_(nullptr),
    useShifted_(realm_.get_shifted_grad_op("pressure")),
    buoyancyWeight_(realm.solutionOptions_->buoyancyPressureStab_ ? 1.0 : 0.0),
    overrideFacePressure_(overrideFacePressure)
{
  // save off fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = metaData.get_field<double>(stk::topology::NODE_RANK, "density");
  interfaceCurvature_ = metaData.get_field<double>(stk::topology::NODE_RANK, "interface_curvature");
  surfaceTension_ = metaData.get_field<double>(stk::topology::NODE_RANK, "surface_tension");
  vof_ = metaData.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  bcPressure_ = metaData.get_field<double>(stk::topology::NODE_RANK, bcPressureName);
  areaWeight_ = metaData.get_field<double>(stk::topology::NODE_RANK, "png_area_weight");
  gravity_ = realm_.solutionOptions_->gravity_;

  NaluEnv::self().naluOutputP0() << "AssembleNodalGradPAWBoundaryAlgorithm Active" << std::endl;
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradPAWBoundaryAlgorithm::execute()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  const int nDim = metaData.spatial_dimension();

  // extract fields
  GenericFieldType *exposedAreaVec = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");

  // ip values; fixed size
  std::vector<double> dpdxBip(nDim);
  std::vector<double> dvofdxBip(nDim);
  double *p_dpdxBip = &dpdxBip[0];
  double *p_dvofdxBip = &dvofdxBip[0];

  // nodal fields to gather; gather everything other than what we are assembling

  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_pressure;
  std::vector<double> ws_face_density;
  std::vector<double> ws_face_kappa;
  std::vector<double> ws_face_sigma;
  std::vector<double> ws_vof;
  std::vector<double> ws_bcPressure;

  // master element
  std::vector<double> ws_face_shape_function;
  std::vector<double> ws_dndx; 
  std::vector<double> ws_det_j;

  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract connected element topology
    b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
    STK_ThrowAssert ( parentTopo.size() == 1 );
    stk::topology theElemTopo = parentTopo[0];

    // volume master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(theElemTopo);
    const int nodesPerElement = meSCS->nodesPerElement_;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_pressure.resize(nodesPerElement);
    ws_face_density.resize(nodesPerFace);
    ws_face_kappa.resize(nodesPerFace);
    ws_face_sigma.resize(nodesPerFace);
    ws_vof.resize(nodesPerElement);
    ws_bcPressure.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);    

    // pointers
    double *p_coordinates = ws_coordinates.data();
    double *p_pressure = ws_pressure.data();
    double *p_face_density = ws_face_density.data();
    double *p_face_kappa = ws_face_kappa.data();
    double *p_face_sigma = ws_face_sigma.data();
    double *p_vof = ws_vof.data();
    double *p_bcPressure = ws_bcPressure.data();
    double *p_face_shape_function = ws_face_shape_function.data();
    double *p_dndx = &ws_dndx[0];

    meFC->shape_fcn(&p_face_shape_function[0]);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulkData.begin_nodes(face);
      int num_face_nodes = bulkData.num_nodes(face);
      // sanity check on num nodes
      STK_ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        p_face_density[ni]    = *stk::mesh::field_data(*density_, node);
        p_face_kappa[ni]    = *stk::mesh::field_data(*interfaceCurvature_, node);
        p_face_sigma[ni]    = *stk::mesh::field_data(*surfaceTension_, node);
        p_bcPressure[ni] = *stk::mesh::field_data(*bcPressure_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec, b, k);

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulkData.begin_elements(face);
      STK_ThrowAssert( bulkData.num_elements(face) == 1 );

      // get element; its face ordinal number and populate face_node_ordinals
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = bulkData.begin_element_ordinals(face)[0];
      const int *face_node_ordinals = meSCS->side_node_ordinals(face_ordinal);

      //======================================
      // gather nodal data off of element
      //======================================
      stk::mesh::Entity const * elem_node_rels = bulkData.begin_nodes(element);
      int num_nodes = bulkData.num_nodes(element);
      // sanity check on num nodes
      STK_ThrowAssert( num_nodes == nodesPerElement );
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        // gather scalars
        p_pressure[ni]    = *stk::mesh::field_data(*pressure_, node);
        p_vof[ni]    = *stk::mesh::field_data(*vof_, node);

        // gather vectors
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // override face nodal pressure (bcPressure will be pressure for all but open)
      if ( overrideFacePressure_ ) {
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const int faceNodeNumber = face_node_ordinals[ic];
          p_pressure[faceNodeNumber] = p_bcPressure[ic];
        }
      }

      // compute dndx
      double scs_error = 0.0;
      if ( useShifted_ )
        meSCS->shifted_face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);
      else
        meSCS->face_grad_op(1, face_ordinal, &p_coordinates[0], &p_dndx[0], &ws_det_j[0], &scs_error);

      // start assembly
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        // nearest node
        const int nn = ipNodeMap[ip];

        stk::mesh::Entity nodeNN = face_node_rels[nn];

        // pointer to fields to assemble
        double *gradPNN = stk::mesh::field_data(*dpdx_, nodeNN);
        double *areaWeightNN = stk::mesh::field_data(*areaWeight_, nodeNN);

        // zero ip values
        double rhoBip = 0.0;
        double sigmaKappaBip = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          p_dpdxBip[j] = 0.0;
          p_dvofdxBip[j] = 0.0;
        }

        // interpolate to bip
        const int ipNpf = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_face_shape_function[ipNpf+ic];
          rhoBip += r*p_face_density[ic];
          sigmaKappaBip += r*p_face_sigma[ic]*p_face_kappa[ic];
        }

        // save off index
        const int ipNdim = ip*nDim;

        // form dpdxBip and dvofdxBip
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          const double pIc = p_pressure[ic];
          const double vofIc = p_vof[ic];
          for ( int j = 0; j < nDim; ++j ) {
            const double dxj = p_dndx[offSetDnDx+j];
            p_dpdxBip[j] += dxj*pIc;
            p_dvofdxBip[j] += dxj*vofIc;
          }
        }

        // assemble to nearest node
        for ( int j = 0; j < nDim; ++j ) {
          const double absArea = std::abs(areaVec[ipNdim+j]);
          const double fac = absArea/rhoBip;
          gradPNN[j] += fac*(p_dpdxBip[j]-buoyancyWeight_*rhoBip*gravity_[j] - sigmaKappaBip*dvofdxBip[j]);
          areaWeightNN[j] += absArea;
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
