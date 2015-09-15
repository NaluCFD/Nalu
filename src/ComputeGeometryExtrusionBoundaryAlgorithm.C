/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeGeometryExtrusionBoundaryAlgorithm.h>

#include <Realm.h>
#include <FieldTypeDef.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <stdexcept>
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ComputeGeometryExtrusionBoundaryAlgorithm - extrude mesh
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeGeometryExtrusionBoundaryAlgorithm::ComputeGeometryExtrusionBoundaryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeGeometryExtrusionBoundaryAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // fields
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  VectorFieldType *haloDxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj");
  VectorFieldType *haloAxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_axj");
  VectorFieldType *edgeAreaVec = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
  ScalarFieldType *dualNodalVolume = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  const int nDim = meta_data.spatial_dimension();

  // zero out halo_axj
  stk::mesh::Selector s_all = stk::mesh::selectUnion(partVec_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * hAxj = stk::mesh::field_data(*haloAxj, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      for ( int j = 0; j < nDim; ++j )
        hAxj[k*nDim+j] = 0.0;
    }
  }

  // proceed with extrusion
  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned );  
  
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
	ib != face_buckets.end() ; ++ib ) {
    
    stk::mesh::Bucket & b = **ib;
    
    // extract master element; hard coded for quad or hex; 
    // quad is always true for 2D while for 3D, either hex or wedge apply
    const stk::topology & theElemTopo = (nDim == 2) ? stk::topology::QUAD_4_2D : stk::topology::HEX_8;
    const int num_face_nodes = (nDim == 2) ? 2 : 4;
    std::vector<int> face_node_ordinals(num_face_nodes);
    
    // extract master element for extruded element type
    MasterElement *meSCS = realm_.get_surface_master_element(theElemTopo);
    MasterElement *meSCV = realm_.get_volume_master_element(theElemTopo);
    
    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    
    // mapping between exposed face and extruded element's overlapping face
    const int *faceNodeOnExtrudedElem = meSCS->faceNodeOnExtrudedElem();

    // mapping between exposed face and extruded element's opposing face
    const int *opposingNodeOnExtrudedElem = meSCS->opposingNodeOnExtrudedElem();
  
    // mapping between exposed face scs ips and halo edge
    const int *faceScsIpOnExtrudedElem = meSCS->faceScsIpOnExtrudedElem();
  
    // mapping between exposed face scs ips and exposed face edge
    const int *faceScsIpOnFaceEdges = meSCS->faceScsIpOnFaceEdges();
    
    // alignment of face:edge ordering and scsip area vector
    const double *edgeAlignedArea = meSCS->edgeAlignedArea();

    // define scratch field
    std::vector<double > ws_coordinates(nodesPerElement*nDim);
    std::vector<double > ws_scs_areav(numScsIp*nDim);
    std::vector<double > ws_scv_volume(nodesPerElement);
    
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get face
      stk::mesh::Entity face = b[k];
      
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const* face_elem_rels = bulk_data.begin_elements(face);
      stk::mesh::ConnectivityOrdinal const* face_elem_ords = bulk_data.begin_element_ordinals(face);
      const int num_elements = bulk_data.num_elements(face);
      ThrowRequire( num_elements == 1 );
      stk::mesh::Entity element = face_elem_rels[0];
      const int face_ordinal = face_elem_ords[0];
      theElemTopo.side_node_ordinals(face_ordinal, face_node_ordinals.begin());
      
      // concentrate on loading up the nodal coordinates for the extruded element
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {
	stk::mesh::Entity node = face_node_rels[ni];
	const double * coords = stk::mesh::field_data(*coordinates, node);
	const double * hDxj = stk::mesh::field_data( *haloDxj, node );
	const int faceNode = faceNodeOnExtrudedElem[face_ordinal*num_nodes + ni];
	const int opposingNode = opposingNodeOnExtrudedElem[face_ordinal*num_nodes + ni];
	const int offSetFN = faceNode*nDim;
	const int offSetON = opposingNode*nDim;
	for ( int j=0; j < nDim; ++j ) {
	  // face node
	  ws_coordinates[offSetFN+j] = coords[j];
	  ws_coordinates[offSetON+j] = coords[j] + hDxj[j];
	}
      }
      
      // compute integration point volume
      double scv_error = 0.0;
      meSCV->determinant(1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);
      
      // compute scs integration point areavec
      double scs_error = 0.0;
      meSCS->determinant(1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);
      
      // assemble nodal dual volume and halo edge area vector
      for ( int ni = 0; ni < num_nodes; ++ni ) {
	stk::mesh::Entity node = face_node_rels[ni];
	double * dualVolume = stk::mesh::field_data( *dualNodalVolume, node );
	const int faceNode = faceNodeOnExtrudedElem[face_ordinal*num_nodes + ni];
	*dualVolume += ws_scv_volume[faceNode];
	// area vector for halo edge; 
	// face ordinal 0 for extruded element has all scs area vectors pointing from face to opposing face
	double * hAxj = stk::mesh::field_data( *haloAxj, node );
	const int scsIp = faceScsIpOnExtrudedElem[face_ordinal*num_nodes + ni];
	for ( int j=0; j < nDim; ++j )
	  hAxj[j] += ws_scs_areav[scsIp*nDim+j];
      }
      
      // deal with edges on the exposed face... 
      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(element);
      
      // face edge relations; if this is 2D then the face is a edge and size is unity
      stk::mesh::Entity const* face_edge_rels = bulk_data.begin_edges(face);
      const int num_face_edges = bulk_data.num_edges(face);
      
      int num_edges = (nDim == 3) ? num_face_edges : 1;
      
      for ( int i = 0; i < num_edges; ++i ) {
	
	// get edge
	stk::mesh::Entity edge = (nDim == 3) ? face_edge_rels[i] : face;
	
	// pointer to edge area vector
	double * av = stk::mesh::field_data(*edgeAreaVec, edge );
	
	// get the relations from edge
	stk::mesh::Entity const* edge_node_rels = bulk_data.begin_nodes(edge);
	const int edge_num_nodes = bulk_data.num_nodes(edge);
	// sanity check on num nodes
	if ( edge_num_nodes != 2 ){
	  throw std::runtime_error("num nodes is not 2");
	}
	
	// extract ip for this edge
	const int scsIp = faceScsIpOnFaceEdges[face_ordinal*num_edges + i];
	
	// correct area for edge and scs area vector from extruded element alignment
	const double alignmentFac = edgeAlignedArea[face_ordinal*num_edges + i];
	
	// area vector points from global low to global high id
	const int leftNode = face_node_ordinals[i];
	const size_t iglob_Lelem = bulk_data.identifier(elem_node_rels[leftNode]);
	const size_t iglob_Ledge = bulk_data.identifier(edge_node_rels[0]);
	
	// determine the sign value for area vector; if Left node is the same,
	// then the element and edge relations are aligned
	const double sign = ( iglob_Lelem == iglob_Ledge ) ? 1.0 : -1.0;
	
	// assemble area vector
	for ( int j=0; j < nDim; ++j )
	  av[j] += ws_scs_areav[scsIp*nDim+j]*alignmentFac*sign;
      }
    }
  }

  // safe to parallel assemble haloAxj here since no one else knows of it...
  std::vector<stk::mesh::FieldBase*> sum_fields(1, haloAxj);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  // dual volume and edge area vector handled elsewhere
}

} // namespace nalu
} // namespace Sierra
