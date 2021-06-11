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
  const std::string bcPressureName)
  : Algorithm(realm, part),
    pressure_(pressure),
    dpdx_(dpdx),
    coordinates_(nullptr),
    density_(nullptr),
    bcPressure_(nullptr),
    areaWeight_(nullptr),
    useShifted_(realm_.get_shifted_grad_op("pressure"))
{
  // save off fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  bcPressure_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, bcPressureName);
  areaWeight_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "png_area_weight");

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
  GenericFieldType *exposedAreaVec = metaData.get_field<GenericFieldType>(metaData.side_rank(), "exposed_area_vector");

  // ip values; fixed size
  std::vector<double> dpdxBip(nDim);
  double *p_dpdxBip = &dpdxBip[0];

  // nodal fields to gather; gather everything other than what we are assembling

  // nodal fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_pressure;
  std::vector<double> ws_face_density;
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
    ThrowAssert ( parentTopo.size() == 1 );
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
    ws_bcPressure.resize(nodesPerFace);
    ws_face_shape_function.resize(numScsBip*nodesPerFace);
    ws_dndx.resize(nDim*numScsBip*nodesPerElement);
    ws_det_j.resize(numScsBip);    

    // pointers
    double *p_coordinates = ws_coordinates.data();
    double *p_pressure = ws_pressure.data();
    double *p_face_density = ws_face_density.data();
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
      ThrowAssert( num_face_nodes == nodesPerFace );
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        // gather scalars
        p_face_density[ni]    = *stk::mesh::field_data(*density_, node);
        p_bcPressure[ni] = *stk::mesh::field_data(*bcPressure_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec, b, k);

      // extract the connected element to this exposed face; should be single in size!
      const stk::mesh::Entity* face_elem_rels = bulkData.begin_elements(face);
      ThrowAssert( bulkData.num_elements(face) == 1 );

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
      ThrowAssert( num_nodes == nodesPerElement );
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        // gather scalars
        p_pressure[ni]    = *stk::mesh::field_data(*pressure_, node);

        // gather vectors
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }

      // override face nodal pressure (bcPressure will be pressure for all but open)
      const bool doIt = false;
      if ( doIt ) {
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
        for ( int j = 0; j < nDim; ++j ) {
          p_dpdxBip[j] = 0.0;
        }

        // interpolate to bip
        const int ipNpf = ip*nodesPerFace;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          rhoBip += p_face_shape_function[ipNpf+ic]*p_face_density[ic];
        }

        // form dpdxBip
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          const double pIc = p_pressure[ic];
          for ( int j = 0; j < nDim; ++j ) {
            p_dpdxBip[j] += p_dndx[offSetDnDx+j]*pIc;
          }
        }

        // assemble to nearest node
        const int ipNdim = ip*nDim;
        for ( int j = 0; j < nDim; ++j ) {
          const double absArea = std::abs(areaVec[ipNdim+j]);
          double fac = absArea/rhoBip;
          gradPNN[j] += fac*p_dpdxBip[j];
          areaWeightNN[j] += absArea;
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
