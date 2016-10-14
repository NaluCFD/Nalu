/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "SuperElement.h"
#include <NaluEnv.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

// edges and faces
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>

// stk_search
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// c++
#include <algorithm>
#include <vector>
#include <stdexcept>

// mpi; for node id consolidation algorithm
#include <mpi.h>

namespace sierra{
namespace naluUnit{

using nalu::NaluEnv;

//==========================================================================
// Class Definition
//==========================================================================
// SuperElement - unit test for super element
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SuperElement::SuperElement()
  : pOrder_(2),
    activateAura_(false),
    currentTime_(0.0),
    resultsFileIndex_(1),
    nDim_(2),
    metaData_(NULL),
    bulkData_(NULL),
    ioBroker_(NULL),
    nodeField_(NULL),
    coordinates_(NULL),
    originalPartName_("block_1"),
    originalSurfacePartName_("surface_1"),
    superElementPartName_("block_1_se"),
    superElementSurfacePartName_("surface_1_se"),
    promotedNodesPartName_("block_1_se_n"),
    edgePartName_("block_1_edges"),
    facePartName_("block_1_faces"),
    verboseOutput_(false),
    originalPart_(NULL),
    originalSurfacePart_(NULL),
    superElementPart_(NULL),
    superSurfacePart_(NULL),
    promotedNodesPart_(NULL),
    edgePart_(NULL),
    facePart_(NULL),
    numberOfEdges_(0),
    numberOfFaces_(0),
    numberOfElements_(0)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SuperElement::~SuperElement()
{
  delete bulkData_;
  delete metaData_;
  delete ioBroker_;
}

void SuperElement::initialize_mesh(stk::ParallelMachine pm)
{
  // news for mesh constructs
  metaData_ = new stk::mesh::MetaData();

  nDim_ = 2;
  metaData_->initialize(nDim_, stk::mesh::entity_rank_names());

  metaData_->declare_part(originalPartName_, stk::topology::ELEM_RANK);
  metaData_->declare_part(originalSurfacePartName_, metaData_->side_rank());
  metaData_->declare_part("surface_2", metaData_->side_rank());
  metaData_->declare_part("surface_3", metaData_->side_rank());
  metaData_->declare_part("surface_4", metaData_->side_rank());
  VectorFieldType& coordfield = metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field(coordfield, metaData_->universal_part(), nDim_, (double*)nullptr);

  bulkData_ = new stk::mesh::BulkData(*metaData_, pm, activateAura_ ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA);
}

void SuperElement::createThreeElemQuadMesh_1proc(stk::mesh::BulkData& mesh)
{
    ASSERT_EQ(2, mesh.mesh_meta_data().spatial_dimension());

    stk::mesh::Part& quadPart = mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4_2D);
    stk::mesh::Part& block1 = *mesh.mesh_meta_data().get_part(originalPartName_);
    stk::mesh::PartVector parts = {&quadPart, &block1};

    stk::mesh::Part& surface1 = *mesh.mesh_meta_data().get_part(originalSurfacePartName_);
    stk::mesh::Part& surface2 = *mesh.mesh_meta_data().get_part("surface_2");
    stk::mesh::Part& surface3 = *mesh.mesh_meta_data().get_part("surface_3");
    stk::mesh::Part& surface4 = *mesh.mesh_meta_data().get_part("surface_4");
    stk::mesh::PartVector surfaceParts = {&surface1};

    mesh.modification_begin();

    stk::mesh::EntityIdVector nodeIds = { 1, 2, 4, 8 };
    stk::mesh::EntityId elemId = 1;
    stk::mesh::Entity elem1 = stk::mesh::declare_element(mesh, parts, elemId, nodeIds);

    nodeIds = { 8, 4, 5, 7 };
    elemId = 2;
    stk::mesh::Entity elem2 = stk::mesh::declare_element(mesh, parts, elemId, nodeIds);

    nodeIds = { 7, 5, 3, 6 };
    elemId = 3;
    stk::mesh::Entity elem3 = stk::mesh::declare_element(mesh, parts, elemId, nodeIds);

    unsigned sideOrdinal = 2;
    mesh.declare_element_side(elem3, sideOrdinal, surfaceParts);
    sideOrdinal = 0;
    surfaceParts[0] = &surface2;
    mesh.declare_element_side(elem1, sideOrdinal, surfaceParts);

    sideOrdinal = 1;
    surfaceParts[0] = &surface3;
    mesh.declare_element_side(elem1, sideOrdinal, surfaceParts);
    mesh.declare_element_side(elem2, sideOrdinal, surfaceParts);
    mesh.declare_element_side(elem3, sideOrdinal, surfaceParts);

    sideOrdinal = 3;
    surfaceParts[0] = &surface4;
    mesh.declare_element_side(elem1, sideOrdinal, surfaceParts);
    mesh.declare_element_side(elem2, sideOrdinal, surfaceParts);
    mesh.declare_element_side(elem3, sideOrdinal, surfaceParts);

    mesh.modification_end();
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void 
SuperElement::execute() 
{
  NaluEnv::self().naluOutputP0() << "Welcome to the SuperElement unit test" << std::endl;
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "SuperElement Quad4 Unit Tests" << std::endl;
  NaluEnv::self().naluOutputP0() << "-----------------------------" << std::endl;

  stk::ParallelMachine pm = NaluEnv::self().parallel_comm();
  
  if (stk::parallel_machine_size(pm) > 1) {
    NaluEnv::self().naluOutputP0() << "!!! SuperElement unit test only runs on 1 proc, exiting early. !!!" << std::endl;
    return;
  }

  initialize_mesh(pm);

  // check to make sure that we are supporting
  if ( nDim_ > 2 || pOrder_ > 2 )
    throw std::runtime_error("Only 2D P=2 is now supported");
  
  // create the part that holds the super element topo (volume and surface)
  declare_super_part();
  declare_super_part_surface();

  declare_edge_part();
  declare_face_part();
  
  // register the fields
  register_fields();
  register_fields_surface();
  
  createThreeElemQuadMesh_1proc(*bulkData_);

  // create the edges and faces on low order part; tmp part(s) to later delete
  create_edges();
  create_faces();
  
  // extract coordinates
  coordinates_ = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  
  // create the parent id maps
  size_of_edges();
  size_of_faces();
  size_of_elements();

  // create nodes
  create_nodes();
  
  // for edges that have multiple processor owners, consolidate ids
  consolidate_node_ids();
  
  // create the element
  create_elements();
  create_elements_surface();
  
  // delete the edges and faces
  delete_edges();
  delete_faces();
  
  // initialize nodal fields; define selector (locally owned and ghosted)
  initialize_fields();
}

//--------------------------------------------------------------------------
//-------- declare_super_part ----------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::declare_super_part()
{
  // set nodes per element; assume quad or hex
  int nodesPerElem = (pOrder_ + 1)*(pOrder_ + 1);
  if ( nDim_ > 2)
    nodesPerElem *= (pOrder_ + 1);

  // create the super topo; how to assign number of sides, faces, etc?
  stk::topology superElemTopo = stk::create_superelement_topology(nodesPerElem);
  
  // create the part all at once
  superElementPart_ = &metaData_->declare_part_with_topology(superElementPartName_, superElemTopo);
  
  // we want this part to show up in the output mesh
  stk::io::put_io_part_attribute(*superElementPart_);
  
  // save off lower order part
  originalPart_ = metaData_->get_part(originalPartName_);
  
  // declare part for nodes
  promotedNodesPart_ = &metaData_->declare_part(promotedNodesPartName_, stk::topology::NODE_RANK);
}

//--------------------------------------------------------------------------
//-------- declare_super_part_surface --------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::declare_super_part_surface()
{
  // now deal with surface
  int nodesPerFace = pOrder_+1;
  if ( nDim_ > 2 )
    nodesPerFace *= (pOrder_+1);

  // create the super topo; how to assign number of sides, faces, etc?
  stk::topology superElemSurfaceTopo = (nDim_ > 2) 
    ? stk::create_superface_topology(nodesPerFace) 
    : stk::create_superedge_topology(nodesPerFace); 
  
  // declare part with superTopo
  superSurfacePart_ = &metaData_->declare_part_with_topology(superElementSurfacePartName_,   superElemSurfaceTopo);
    
  // we want this part to show up in the output mesh... may not work as expected
  const bool outputIt = false;
  if ( outputIt )
    stk::io::put_io_part_attribute(*superSurfacePart_);
    
  // save off lower order part
  originalSurfacePart_ = metaData_->get_part(originalSurfacePartName_);
}

//--------------------------------------------------------------------------
//-------- declare_edge_part -----------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::declare_edge_part()
{
  edgePart_ = &metaData_->declare_part(edgePartName_, stk::topology::NODE_RANK);
}

//--------------------------------------------------------------------------
//-------- declare_face_part -----------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::declare_face_part()
{
  if ( nDim_ == 3 )
    facePart_ = &metaData_->declare_part(edgePartName_, stk::topology::FACE_RANK);
}

//--------------------------------------------------------------------------
//-------- create_edges ----------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::create_edges()
{
  stk::mesh::create_edges(*bulkData_, stk::mesh::Selector(*originalPart_), edgePart_);
}

//--------------------------------------------------------------------------
//-------- delete_edges ----------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::delete_edges()
{
  // extract selected edges; avoid bucket loop since it gets out of date
  stk::mesh::BucketVector const& edge_buckets = bulkData_->get_buckets( stk::topology::EDGE_RANK,  *edgePart_);
  std::vector<stk::mesh::Entity> edges;
  stk::mesh::get_selected_entities( *edgePart_ , edge_buckets, edges);

  // delete upward relations
  bulkData_->modification_begin();
  for (unsigned ii=0; ii < edges.size(); ++ii) {

    if (!bulkData_->is_valid(edges[ii]))
      throw std::runtime_error("bad edge 1");

    // upward edge->element
    unsigned num_elems = bulkData_->num_elements(edges[ii]);
    if ( num_elems > 0 ) {
      stk::mesh::Entity const * const edge_elems = bulkData_->begin_elements(edges[ii]);
      stk::mesh::ConnectivityOrdinal const* edge_elem_ordinals = bulkData_->begin_element_ordinals(edges[ii]);
      for ( unsigned k = 0; k < num_elems; ++k ) {
        stk::mesh::Entity to_rel = edge_elems[k];
        stk::mesh::RelationIdentifier to_id = edge_elem_ordinals[k];
        bool del = bulkData_->destroy_relation( to_rel, edges[ii], to_id);
        if (!del)
          throw std::runtime_error("delete_edges failed to delete up relation (element)");
      }
    }
    
    // upward edge->face
    unsigned num_faces = bulkData_->num_faces(edges[ii]);
    if ( num_faces > 0 ) {
      stk::mesh::Entity const * const edge_faces = bulkData_->begin_faces(edges[ii]);
      stk::mesh::ConnectivityOrdinal const* edge_face_ordinals = bulkData_->begin_face_ordinals(edges[ii]);
      for ( unsigned k = 0; k < num_faces; ++k ) {
        stk::mesh::Entity to_rel = edge_faces[k];
        stk::mesh::RelationIdentifier to_id = edge_face_ordinals[k];
        bool del = bulkData_->destroy_relation( to_rel, edges[ii], to_id);
        if (!del)
          throw std::runtime_error("delete_edges failed to delete up relation (face)");
      }
    }

    // good to destroy
    if (bulkData_->is_valid(edges[ii]) && bulkData_->bucket(edges[ii]).owned()) {
      const bool successfullyDestroyed = bulkData_->destroy_entity( edges[ii] );   
      if ( !successfullyDestroyed )
        throw std::runtime_error("destroy did not happen");
    }
  }  
  bulkData_->modification_end();
  
  // check now...
  size_t numberOfEdgesInOriginalPart = 0;
  stk::mesh::Selector s_edge_orig_part = stk::mesh::Selector(*originalPart_);  
  stk::mesh::BucketVector const& edge_buckets_orig_part =
    bulkData_->get_buckets(stk::topology::EDGE_RANK, s_edge_orig_part );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets_orig_part.begin();
        ib != edge_buckets_orig_part.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();    
    // increment size
    numberOfEdgesInOriginalPart += length; 
  }

  // second, check the size of edges on the edge part
  size_t numberOfEdgesInEdgePart = 0;
  stk::mesh::Selector s_edge_edge_part = stk::mesh::Selector(*edgePart_);
  stk::mesh::BucketVector const& edge_buckets_edge_part =
    bulkData_->get_buckets(stk::topology::EDGE_RANK, s_edge_edge_part );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets_edge_part.begin();
        ib != edge_buckets_edge_part.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();    
    // increment size
    numberOfEdgesInEdgePart += length; 
  }
  
  // parallel sum
  size_t l_sum[2] = {numberOfEdgesInOriginalPart, numberOfEdgesInEdgePart};
  size_t g_sum[2] = {};
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_sum(comm, l_sum, g_sum, 2);
  // reassign
  numberOfEdgesInOriginalPart = g_sum[0];
  numberOfEdgesInEdgePart = g_sum[1];

  // there should be two removed from the edge part (zero left); 
  // there should be eight remaining in the original part
  const bool testEdgesInOrigPart = numberOfEdgesInOriginalPart == 8 ? true : false;
  const bool testEdgesInEdgePart = numberOfEdgesInEdgePart == 0 ? true : false;
    
  if ( testEdgesInOrigPart )
    NaluEnv::self().naluOutputP0() << "Remaining Edges(OP) Test    PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Remaining Edges(OP) Test    FAILED" << std::endl;
  
  if ( testEdgesInEdgePart )
    NaluEnv::self().naluOutputP0() << "Remaining Edges(EP) Test    PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Remaining Edges(EP) Test    FAILED" << std::endl;
}

//--------------------------------------------------------------------------
//-------- create_faces ------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::create_faces()
{
  if ( nDim_ == 3 )
    stk::mesh::create_faces(*bulkData_, stk::mesh::Selector(*originalPart_), facePart_);
}

//--------------------------------------------------------------------------
//-------- delete_faces ------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::delete_faces()
{
  // delete it; not yet
}

//--------------------------------------------------------------------------
//-------- size_of_edges ---------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::size_of_edges()
{ 
  // size edge count
  numberOfEdges_ = 0;

  // selector based on locally owned and shared edges
  stk::mesh::Selector s_edge = stk::mesh::Selector(*originalPart_);

  stk::mesh::BucketVector const& all_edge_buckets =
    bulkData_->buckets(stk::topology::EDGE_RANK);
  stk::mesh::BucketVector const& edge_buckets =
    bulkData_->get_buckets(stk::topology::EDGE_RANK, s_edge );
EXPECT_EQ(all_edge_buckets.size(), edge_buckets.size());
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();    

    // increment size
    numberOfEdges_ += length;
      
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {           
      // finally, determine the the set of processors that this edge touches; defines commmunication pattern for nodes
      std::vector<int> sharedProcsEdge;
      bulkData_->comm_shared_procs({stk::topology::EDGE_RANK, bulkData_->identifier(b[k])}, sharedProcsEdge);
      std::sort(sharedProcsEdge.begin(), sharedProcsEdge.end());
      sharedProcsEdge_.push_back(sharedProcsEdge);
    }
  }

  if (verboseOutput_ )
    NaluEnv::self().naluOutputP0() << "size of edges: " << numberOfEdges_ << std::endl;

  EXPECT_EQ(10, numberOfEdges_);
  const bool testEdge = numberOfEdges_ == 10 ? true : false;
  if ( testEdge )
    NaluEnv::self().naluOutputP0() << "Total Edge Count Test       PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Total Edge Count Test       FAILED" << std::endl;
}

//--------------------------------------------------------------------------
//-------- size_of_faces ---------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::size_of_faces()
{
  // size number of faces; not ready for 3D....
  numberOfFaces_ = 0;

  // define some common selectors
  stk::mesh::Selector s_face = stk::mesh::Selector(*originalPart_);

  stk::mesh::BucketVector const& face_buckets =
    bulkData_->get_buckets(stk::topology::FACE_RANK, s_face );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();    
    // increment size
    numberOfFaces_ += length;
  }
  
  if (verboseOutput_ )
    NaluEnv::self().naluOutputP0() << "size of faces: " << numberOfFaces_ << std::endl;

  if ( numberOfFaces_ > 0 )
    throw std::runtime_error("size_of_faces: greater than zero; 3D not ready for prime time");  
}

//--------------------------------------------------------------------------
//-------- size_of_elements ------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::size_of_elements()
{
  // size element count
  numberOfElements_ = 0;

  // define some common selectors; want locally owned here
  stk::mesh::Selector s_elem = metaData_->locally_owned_part()
    & stk::mesh::Selector(*originalPart_);

  stk::mesh::BucketVector const& elem_buckets =
    bulkData_->get_buckets(stk::topology::ELEMENT_RANK, s_elem );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    // increment size
    numberOfElements_ += length;   
  }

  if (verboseOutput_ )
    NaluEnv::self().naluOutputP0() << "size of elems: " << numberOfElements_ << std::endl;

  const bool testElem = numberOfElements_ == 3 ? true : false;
  if ( testElem )
    NaluEnv::self().naluOutputP0() << "Total Elem Count Test       PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Total Elem Count Test       FAILED" << std::endl;
}

//--------------------------------------------------------------------------
//-------- create_nodes ----------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::create_nodes()
{
  // count the number of promoted nodal ids required; based on parentIds size (which has been sorted)
  const int pM1Order = pOrder_ - 1;
  const int pElemFac = std::pow(pM1Order, nDim_);
  const int pEdgeFac = pM1Order;
  const int pFaceFac = pM1Order*pM1Order;
  
  const int numPromotedNodes
    = numberOfEdges_*pEdgeFac
    + numberOfFaces_*pFaceFac 
    + numberOfElements_*pElemFac;

  // okay, now ask
  bulkData_->modification_begin();

  // generate new ids; number of points is simple for now... all of the extra nodes from P=1 to P=2
  stk::mesh::EntityIdVector availableNodeIds(numPromotedNodes);
  bulkData_->generate_new_ids(stk::topology::NODE_RANK, numPromotedNodes, availableNodeIds);
  
  // declare the entity on this rank (rank is determined by calling declare_entity on this rank)
  for (int i = 0; i < numPromotedNodes; ++i) {
    stk::mesh::Entity theNode 
      = bulkData_->declare_entity(stk::topology::NODE_RANK, availableNodeIds[i], *promotedNodesPart_);
    promotedNodesVec_.push_back(theNode);
  }

  bulkData_->modification_end();
      
  // fill in std::map<stk::mesh::EntityIdVector, stk::mesh::Entity > parentNodesMap_
  int promotedNodesVecCount = 0;

  // edge selectors; locally owned and shared edges
  stk::mesh::Selector s_edge = stk::mesh::Selector(*originalPart_);

  stk::mesh::BucketVector const& edge_buckets =
    bulkData_->get_buckets(stk::topology::EDGE_RANK, s_edge );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {           
      
      // get edge and edge id
      stk::mesh::Entity edge = b[k];
      const stk::mesh::EntityId edgeId = bulkData_->identifier(edge);

      // extract node relations from the edge and node count
      stk::mesh::Entity const * edge_node_rels =  bulkData_->begin_nodes(edge);
      const int nodesPerEdge = b.num_nodes(k);

      // sanity check on number or nodes
      ThrowAssert( 2 == nodesPerEdge );

      // left, right and center node along the edge
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];
      stk::mesh::Entity nodeC = promotedNodesVec_[promotedNodesVecCount];
      
      const double * edgeCoordsL = stk::mesh::field_data(*coordinates_, nodeL);
      const double * edgeCoordsR = stk::mesh::field_data(*coordinates_, nodeR);
      double * edgeCoordsC = stk::mesh::field_data(*coordinates_, nodeC);
      
      // find mean distance between the nodes
      for ( int j = 0; j < nDim_; ++j )
        edgeCoordsC[j] = (edgeCoordsL[j] + edgeCoordsR[j])*0.5;
      
      // store off map; for now, only one node per edge (P=2)
      std::vector<stk::mesh::Entity> edgeNodesVec;
      edgeNodesVec.push_back(nodeC);
      parentEdgeNodesMap_[edgeId] = edgeNodesVec;
      
      ++promotedNodesVecCount;
    }
  }

  // fill in faces
  stk::mesh::Selector s_face = stk::mesh::Selector(*originalPart_);

  stk::mesh::BucketVector const& face_buckets =
    bulkData_->get_buckets(stk::topology::FACE_RANK, s_face );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {           
      
      // get face and face id
      stk::mesh::Entity face = b[k];
      const stk::mesh::EntityId faceId = bulkData_->identifier(face);

      // extract node relations from the face and node count
      stk::mesh::Entity const * face_node_rels =  bulkData_->begin_nodes(face);
      const int numNodes = b.num_nodes(k);

      // extract element center node; FIXME: P=2
      stk::mesh::Entity nodeC = promotedNodesVec_[promotedNodesVecCount];

      double * elemCoordsC = stk::mesh::field_data(*coordinates_, nodeC);
    
      // hacked center coords calulation
      std::vector<double>tmpCoord(nDim_,0.0);
      for ( int ni = 0; ni < numNodes; ++ni ) {
        stk::mesh::Entity theNode = face_node_rels[ni];
        double * elemNodeCoords = stk::mesh::field_data(*coordinates_, theNode);
        for ( int i = 0; i < nDim_; ++i )
          tmpCoord[i] += elemNodeCoords[i]/(double)numNodes;
      }
      
      for ( int i = 0; i < nDim_; ++i )
        elemCoordsC[i] = tmpCoord[i];

      // FIXME: P=2
      std::vector<stk::mesh::Entity> faceNodesVec;
      faceNodesVec.push_back(nodeC);      
      parentFaceNodesMap_[faceId] = faceNodesVec;
      
      ++promotedNodesVecCount;
    }
  }

  // element selectors; locally owned only
  stk::mesh::Selector s_elem = metaData_->locally_owned_part()
    & stk::mesh::Selector(*originalPart_);
  
  stk::mesh::BucketVector const& elem_buckets =
    bulkData_->get_buckets(stk::topology::ELEMENT_RANK, s_elem );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get element and element id
      stk::mesh::Entity elem = b[k];
      const stk::mesh::EntityId elemId = bulkData_->identifier(elem);
      
      // extract node relations
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);

      // iterate over the nodes in the element
      int numNodes = b.num_nodes(k);

      // extract element center noce
      stk::mesh::Entity nodeC = promotedNodesVec_[promotedNodesVecCount];
      double * elemCoordsC = stk::mesh::field_data(*coordinates_, nodeC);
    
      // hacked center coords calulation
      std::vector<double>tmpCoord(nDim_,0.0);
      for ( int ni = 0; ni < numNodes; ++ni ) {
        stk::mesh::Entity theNode = node_rels[ni];
        double * elemNodeCoords = stk::mesh::field_data(*coordinates_, theNode);
        for ( int i = 0; i < nDim_; ++i )
          tmpCoord[i] += elemNodeCoords[i]/(double)numNodes;
      }
      
      for ( int i = 0; i < nDim_; ++i )
        elemCoordsC[i] = tmpCoord[i];

      std::vector<stk::mesh::Entity> elemNodesVec;
      elemNodesVec.push_back(nodeC);
      parentElemNodesMap_[elemId] = elemNodesVec;
      
      ++promotedNodesVecCount;
    }  
  }
}

//--------------------------------------------------------------------------
//-------- consolidate_node_ids --------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::consolidate_node_ids()
{
  // generic iterator for parentNodesMap_
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> >::iterator iterFindMap;

  // edge selectors; shared edges only in this part
  stk::mesh::Selector s_edge = stk::mesh::Selector(*originalPart_) & metaData_->globally_shared_part();
  
  stk::mesh::BucketVector const& edge_buckets =
    bulkData_->get_buckets(stk::topology::EDGE_RANK, s_edge );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {           
      
      // get edge and edge id
      stk::mesh::Entity edge = b[k];
      const stk::mesh::EntityId edgeId = bulkData_->identifier(edge);

      // find the super node off the edge
      std::vector<stk::mesh::Entity> edgeNodesVec;
      iterFindMap = parentEdgeNodesMap_.find(edgeId);
      if ( iterFindMap != parentEdgeNodesMap_.end() ) {
        edgeNodesVec = iterFindMap->second;
      }
      else {
        throw std::runtime_error("Could not find the node(s) beloging to low order edge id: " + edgeId );
      }
     
      // determine how many procs touch this edge
      stk::mesh::EntityKey key = bulkData_->entity_key(edge);
      std::vector<int> procs; // should be at least one in size..
      bulkData_->comm_procs(key, procs);

      // create the struct, fill in the data and push back
      for ( size_t iProc = 0; iProc < procs.size(); ++iProc ) {
        for ( size_t iNode = 0; iNode < edgeNodesVec.size(); ++iNode ) {
          EntityNodeSharing enShare;
          enShare.edgeId_ = edgeId;
          enShare.globalNodeId_ = bulkData_->identifier(edgeNodesVec[iNode]);
          enShare.localIndex_ = iNode;
          edgeNodeSharingMap_[procs[iProc]].push_back(enShare);
        }
      }
    }
  }
   
  // setup
  int numNeighbors = edgeNodeSharingMap_.size();
  std::vector<MPI_Request> sendRequests(numNeighbors);
  std::vector<MPI_Request> recvRequests(numNeighbors);

  // find the size of a struct
  EntityNodeSharing enScratch;
  const int sizeofENS = sizeof(enScratch);
  
  int counter = 0;
  std::map<int, std::vector<EntityNodeSharing> >::const_iterator iterM;
  for ( iterM = edgeNodeSharingMap_.begin(); iterM != edgeNodeSharingMap_.end(); ++iterM ) {
    const int senderProc = iterM->first;
    const std::vector<EntityNodeSharing> sendVec = iterM->second;
    const int buffSize  = sendVec.size()*sizeofENS;
    edgeNodeSharingOffProcMap_[senderProc].resize(sendVec.size());
    MPI_Irecv(edgeNodeSharingOffProcMap_[senderProc].data(), buffSize, MPI_BYTE, senderProc,
              MPI_ANY_TAG, NaluEnv::self().parallel_comm(), &recvRequests[counter]);
    counter++;
  }
  
  counter = 0;
  for ( iterM = edgeNodeSharingMap_.begin(); iterM != edgeNodeSharingMap_.end(); ++iterM ) {
    const int destProc = iterM->first;
    const std::vector<EntityNodeSharing> sendVec = iterM->second;
    const int buffSize  = sendVec.size()*sizeofENS;
    const int whatTagIsThis = 0;
    MPI_Isend(sendVec.data(), buffSize, MPI_BYTE, destProc, whatTagIsThis, NaluEnv::self().parallel_comm(), &sendRequests[counter]);
    counter++;
  }
  
  std::vector<MPI_Status> sendStatus(numNeighbors);
  std::vector<MPI_Status> recvStatus(numNeighbors);
     
  MPI_Waitall(numNeighbors, recvRequests.data(), recvStatus.data() );
  MPI_Waitall(numNeighbors, sendRequests.data(), sendStatus.data() ); 
  
  bulkData_->modification_begin();

  for ( iterM = edgeNodeSharingOffProcMap_.begin(); iterM != edgeNodeSharingOffProcMap_.end(); ++iterM ) {
    const std::vector<EntityNodeSharing> receiveBuffer = iterM->second;
    for ( size_t k = 0; k < receiveBuffer.size(); ++k) {
      stk::mesh::EntityId theId = receiveBuffer[k].globalNodeId_;
      stk::mesh::EntityId edgeId = receiveBuffer[k].edgeId_;
      int localIndex = receiveBuffer[k].localIndex_;
      iterFindMap = parentEdgeNodesMap_.find(edgeId);
      if ( iterFindMap != parentEdgeNodesMap_.end() ) {
        stk::mesh::Entity foundNode = iterFindMap->second[localIndex];
        stk::mesh::EntityId foundNodeId = bulkData_->identifier(foundNode);
        stk::mesh::EntityId minId = std::min(foundNodeId, theId);
        if ( foundNodeId != minId ) {
          bulkData_->change_entity_id(minId,foundNode);
        }
      }
      else {
        throw std::runtime_error("Could not find the node(s) beloging to low order edge id: " + edgeId );
      } 
    }
  }

  bulkData_->modification_end();

 
  // check
  const bool parallelCheck = false;
  const bool sendItOut = false;
  if ( parallelCheck ) {

    stk::mesh::Selector s_all_edge = stk::mesh::Selector(*originalPart_);
  
    stk::mesh::BucketVector const& all_edge_buckets =
      bulkData_->get_buckets(stk::topology::EDGE_RANK, s_all_edge );
    for ( stk::mesh::BucketVector::const_iterator ib = all_edge_buckets.begin();
          ib != all_edge_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {           
        
        // get edge and edge id
        stk::mesh::Entity edge = b[k];
        const stk::mesh::EntityId edgeId = bulkData_->identifier(edge);
        
        // find the super node off the edge
        std::vector<stk::mesh::Entity> edgeNodesVec;
        iterFindMap = parentEdgeNodesMap_.find(edgeId);
        if ( iterFindMap != parentEdgeNodesMap_.end() ) {
          edgeNodesVec = iterFindMap->second;
        }
        else {
          throw std::runtime_error("Could not find the node(s) beloging to low order edge id: " + edgeId );
        }
        
        if ( sendItOut )
          NaluEnv::self().naluOutput() << "edge Id: " << edgeId << std::endl;
        for ( size_t iNode = 0; iNode < edgeNodesVec.size(); ++iNode) {
          stk::mesh::EntityId nodeId = bulkData_->identifier(edgeNodesVec[iNode]);
          if ( sendItOut )
            NaluEnv::self().naluOutput() << "node id: " << nodeId << std::endl;
        } 
      }
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- create_elements -------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::create_elements()
{
  // define some common selectors; want locally owned here
  stk::mesh::Selector s_elem = metaData_->locally_owned_part()
    & stk::mesh::Selector(*originalPart_);
  
  stk::mesh::BucketVector const& elem_buckets =
    bulkData_->get_buckets(stk::topology::ELEMENT_RANK, s_elem );

  // elements and assign the new node connectivity
  bulkData_->modification_begin();
  
  // generate new ids; one per bucket loop
  stk::mesh::EntityIdVector availableElemIds(numberOfElements_);
  bulkData_->generate_new_ids(stk::topology::ELEM_RANK, numberOfElements_, availableElemIds);

  // generic iterator for parentNodesMap_
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> >::iterator iterFindMap;
  
  // declare id counter
  size_t availableElemIdCounter = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
         ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // define the vector that will hold the connected nodes for this element
      stk::mesh::EntityIdVector connectedNodesIdVec;
      
      // get element and element id
      stk::mesh::Entity elem = b[k];
      const stk::mesh::EntityId elemId = bulkData_->identifier(elem);
      
      // extract node relations amd mpde count
      stk::mesh::Entity const * elem_node_rels =  bulkData_->begin_nodes(elem);
      int numElemNodes = b.num_nodes(k);

      // first, standard nodes on the lower order elements
      for ( int ni = 0; ni < numElemNodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        connectedNodesIdVec.push_back(bulkData_->identifier(node));
      }
      
      // second, nodes in the center of the edges
      stk::mesh::Entity const * elem_edge_rels = bulkData_->begin_edges(elem);
      int numEdges = b.num_edges(k);
 
      for ( int ne = 0; ne < numEdges; ++ne ) {
         
        // extract the edge and edge id
        stk::mesh::Entity edge = elem_edge_rels[ne];
        stk::mesh::EntityId edgeId = bulkData_->identifier(edge);
        
        // find the edge centroid node(s)
        iterFindMap = parentEdgeNodesMap_.find(edgeId);
        std::vector<stk::mesh::Entity> edgeNodesVec;
        if ( iterFindMap != parentEdgeNodesMap_.end() ) {
          edgeNodesVec = iterFindMap->second;
        }
        else {
          throw std::runtime_error("Could not find the node beloging to edge vector");
        }
        for ( size_t iNode = 0; iNode < edgeNodesVec.size(); ++iNode )
          connectedNodesIdVec.push_back(bulkData_->identifier(edgeNodesVec[iNode]));
      }
                  
      // last, nodes in the center of the element; find the element centroid node(s)
      iterFindMap = parentElemNodesMap_.find(elemId);
      std::vector<stk::mesh::Entity> elemNodesVec;
      if ( iterFindMap != parentElemNodesMap_.end() ) {
        elemNodesVec = iterFindMap->second;
      }
      else {
        throw std::runtime_error("Could not find the node beloging to element vector");
      }
      
      // FIXME: assumes P=1
      for (size_t iNode = 0; iNode < elemNodesVec.size(); ++iNode)
        connectedNodesIdVec.push_back(bulkData_->identifier(elemNodesVec[iNode]));
      
      // all done with element, edge and face node connectivitoes; create the element
      stk::mesh::Entity theElem
        = stk::mesh::declare_element(*bulkData_, *superElementPart_,
                                     availableElemIds[availableElemIdCounter],
                                     connectedNodesIdVec);
      
      // push back to map
      superElementElemMap_[elemId] = theElem;
      availableElemIdCounter++;
    }
  }
  
  bulkData_->modification_end();
}

//--------------------------------------------------------------------------
//-------- create_elements_surface ----------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::create_elements_surface()
{
  // find total number of locally owned elements
  size_t numNewSurfaceElem = 0;

  // generic iterator for parentNodesMap_ and placeholder for the found element
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> >::iterator iterFindMap;
  stk::mesh::Entity foundElem;
  
  // define vector of parent topos; should always be UNITY in size
  std::vector<stk::topology> parentTopo;

  // define some common selectors
  stk::mesh::Selector s_low_order = metaData_->locally_owned_part()
  & stk::mesh::Selector(*originalSurfacePart_);
  
  stk::mesh::BucketVector const& face_buckets =
    bulkData_->get_buckets(metaData_->side_rank(), s_low_order );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    numNewSurfaceElem += length;
  }
  
  // now loop over elements and assign the new node connectivity
  bulkData_->modification_begin();
  
  // generate new ids; one per bucket loop
  stk::mesh::EntityIdVector availableSurfaceElemIds(numNewSurfaceElem);
  bulkData_->generate_new_ids(metaData_->side_rank(), numNewSurfaceElem, availableSurfaceElemIds);
  
  // declare id counter
  size_t availableSurfaceElemIdCounter = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    // extract buckets face and element topology; not sure if we need this yet
    /*
      stk::topology thisBucketsTopo = b.topology();
      b.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
      ThrowAssert ( parentTopo.size() == 1 );
      stk::topology theElemTopo = parentTopo[0];
    */

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get face and its id
      stk::mesh::Entity face = b[k];
      const stk::mesh::EntityId faceId = bulkData_->identifier(face);
    
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = bulkData_->begin_elements(face);
      ThrowAssert( bulkData_->num_elements(face) == 1 );
      
      // get element, element id, and face ordinal number
      stk::mesh::Entity elem = face_elem_rels[0];
      const stk::mesh::EntityId elemId = bulkData_->identifier(elem);
      const int faceOrdinal = bulkData_->begin_element_ordinals(face)[0];
      
      // set local node count used for face:node relations
      int localNodeCount = 0;

      // declare the super face (with side rank) with current surface id
      stk::mesh::Entity superFace 
        = bulkData_->declare_entity(metaData_->side_rank(), 
                                    availableSurfaceElemIds[availableSurfaceElemIdCounter], 
                                    *superSurfacePart_);
      
      // find the super element
      std::map<stk::mesh::EntityId, stk::mesh::Entity>::iterator iterFindElem;
      iterFindElem = superElementElemMap_.find(elemId);
      if ( iterFindElem != superElementElemMap_.end() ) {
        foundElem = iterFindElem->second;
      }
      else {
        throw std::runtime_error("Could not find the super element beloging to low order element id: " + elemId );
      }
      
      // first, declare face:element relation
      bulkData_->declare_relation(foundElem, superFace, faceOrdinal);

      // next, face:node relations; parent element's face node first
      stk::mesh::Entity const * face_node_rels = bulkData_->begin_nodes(face);
      int num_face_nodes = bulkData_->num_nodes(face);
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        bulkData_->declare_relation(superFace, node, localNodeCount);
        localNodeCount++;
      }
      
      // find the connected edges to this face
      stk::mesh::Entity const* face_edge_rels = bulkData_->begin_edges(face);

      const int num_face_edges = bulkData_->num_edges(face);
      
      int num_edges = (nDim_ == 3) ? num_face_edges : 1;
      
      for ( int i = 0; i < num_edges; ++i ) {
	
	// get edgeId
	const stk::mesh::EntityId edgeId = (nDim_ == 3) ? bulkData_->identifier(face_edge_rels[i]) : faceId;

        // find the super node off the edge
        std::vector<stk::mesh::Entity> edgeNodesVec;
        iterFindMap = parentEdgeNodesMap_.find(edgeId);
        if ( iterFindMap != parentEdgeNodesMap_.end() ) {
          edgeNodesVec = iterFindMap->second;
        }
        else {
          throw std::runtime_error("Could not find the node(s) beloging to low order edge id: " + edgeId );
        }
        
        // add nodes on each edge relation
        for ( size_t j = 0; j < edgeNodesVec.size(); ++j) {
          bulkData_->declare_relation(superFace, edgeNodesVec[j], localNodeCount);
          localNodeCount++;
        }
      }

      // now faces
      if ( nDim_ == 3 ) {  
        std::vector<stk::mesh::Entity> faceNodesVec;
        // find the super node off the face
        iterFindMap = parentFaceNodesMap_.find(faceId);
        if ( iterFindMap != parentFaceNodesMap_.end() ) {
          faceNodesVec = iterFindMap->second;
        }
        else {
          throw std::runtime_error("Could not find the node(s) beloging to low order edge id: " + faceId );
        }     
      
        // add final set of relations
        for ( size_t j = 0; j < faceNodesVec.size(); ++j) {
          bulkData_->declare_relation(superFace, faceNodesVec[j], localNodeCount);
          localNodeCount++;
        }
      }

      // increment available ids
      availableSurfaceElemIdCounter++;
    }
  }
  
  bulkData_->modification_end();

  //=========================================
  // now check surface 1
  //=========================================

  // set some gold standards
  std::vector<stk::mesh::EntityId > nodeIdGold(3);
  std::vector<stk::mesh::EntityId > nodeIdCheck;
  nodeIdGold[0] = 3;
  nodeIdGold[1] = 6;
  nodeIdGold[2] = 9;
  const int faceIdGold = 37;
  const int nodesPerFaceGold = 3;
  const int faceOrdinalGold = 2;
  const int elemIdGold = 6;
  
  // logics
  bool nodesPerFaceCheck = true;
  bool faceOrdinalCheck = true;
  bool faceIdCheck = true;
  bool elemIdCheck = true; 
  bool testNodeSizeCheck = true;
  bool testNodeIdCheck = true;

  // define some common selectors
  stk::mesh::Selector s_high_order = metaData_->locally_owned_part()
  & stk::mesh::Selector(*superSurfacePart_);
  
  stk::mesh::BucketVector const& high_face_buckets =
    bulkData_->get_buckets(metaData_->side_rank(), s_high_order );
  for ( stk::mesh::BucketVector::const_iterator ib = high_face_buckets.begin();
        ib != high_face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
        
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get face and its id
      stk::mesh::Entity face = b[k];
      const stk::mesh::EntityId faceId = bulkData_->identifier(face);
      const int nodesPerFace = b.topology().num_nodes();
      
      EXPECT_EQ(faceIdGold, faceId);
      if ( faceId != faceIdGold )
        faceIdCheck = false;

      if ( nodesPerFace != nodesPerFaceGold) 
        nodesPerFaceCheck = false;
      
      if ( verboseOutput_ )
        NaluEnv::self().naluOutputP0()  
          << "faceId and nodes per face " << faceId << " " << nodesPerFace << std::endl;
      
      // extract the connected element to this exposed face; should be single in size!
      stk::mesh::Entity const * face_elem_rels = bulkData_->begin_elements(face);
      ThrowAssert( bulkData_->num_elements(face) == 1 );
      
      // get element, element id, and face ordinal number
      stk::mesh::Entity elem = face_elem_rels[0];
      const stk::mesh::EntityId elemId = bulkData_->identifier(elem);
      const int faceOrdinal = bulkData_->begin_element_ordinals(face)[0];

      EXPECT_EQ(faceOrdinalGold, faceOrdinal);
      if ( faceOrdinal != faceOrdinalGold )
        faceOrdinalCheck = false;

      EXPECT_EQ(elemIdGold, elemId);
      if ( elemId != elemIdGold )
        elemIdCheck = false;

      if ( verboseOutput_ ) {
        NaluEnv::self().naluOutputP0() << "elemId, face ordinal and nodesPerFace " << elemId << " " 
                                       << faceOrdinal << " " << nodesPerFace << std::endl;
      }

      stk::mesh::Entity const * face_node_rels = bulkData_->begin_nodes(face);
      int num_face_nodes = bulkData_->num_nodes(face);
      // sanity check on num nodes
      if ( num_face_nodes != nodesPerFace )
        nodesPerFaceCheck = false;
      for ( int ni = 0; ni < num_face_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        const stk::mesh::EntityId nodeId = bulkData_->identifier(node);
        if ( verboseOutput_ ) {
          NaluEnv::self().naluOutputP0() << " surface check node id " << nodeId << std::endl;
        }
        nodeIdCheck.push_back(nodeId);
      }
    }
  }

  if ( nodeIdCheck.size() != nodeIdGold.size()) {
    testNodeSizeCheck = false;
    testNodeIdCheck = false;
  }
  else {
    for ( size_t k = 0; k < nodeIdGold.size(); ++k ) {
      EXPECT_EQ(nodeIdGold[k], nodeIdCheck[k]);
      if ( nodeIdGold[k] != nodeIdCheck[k] ) {
        testNodeIdCheck = false;
      }
    } 
  }

  if ( testNodeSizeCheck )
    NaluEnv::self().naluOutputP0() << "Surface node size Test      PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Surface node size Test      FAILED" << std::endl;

  if ( testNodeIdCheck )
    NaluEnv::self().naluOutputP0() << "Surface node id Test        PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Surface node id Test        FAILED" << std::endl;  

 if ( nodesPerFaceCheck )
    NaluEnv::self().naluOutputP0() << "Surface npf Test            PASSED" << std::endl;
 else
    NaluEnv::self().naluOutputP0() << "Surface npf Test            FAILED" << std::endl;

 if ( faceOrdinalCheck )
    NaluEnv::self().naluOutputP0() << "Surface face ordinal Test   PASSED" << std::endl;
 else
    NaluEnv::self().naluOutputP0() << "Surface face ordinal Test   FAILED" << std::endl;

 if ( faceIdCheck )
    NaluEnv::self().naluOutputP0() << "Surface face Id Test        PASSED" << std::endl;
 else
    NaluEnv::self().naluOutputP0() << "Surface face Id Test        FAILED" << std::endl;

 if ( elemIdCheck )
    NaluEnv::self().naluOutputP0() << "Surface elem Id Test        PASSED" << std::endl;
 else
    NaluEnv::self().naluOutputP0() << "Surface elem Id Test        FAILED" << std::endl;
}
  
//--------------------------------------------------------------------------
//-------- register_fields -------------------------------------------------
//--------------------------------------------------------------------------
void 
SuperElement::register_fields()
{        
  // declare and put nodal field on part
  nodeField_ = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "node_field"));
    
  // put them on the part
  stk::mesh::put_field(*nodeField_, *superElementPart_);
  
  // declare and put element field on part
  elemField_ = &(metaData_->declare_field<GenericFieldType>(stk::topology::ELEM_RANK, "elem_field"));
  stk::mesh::put_field(*elemField_, *superElementPart_, 1);
}

//--------------------------------------------------------------------------
//-------- register_fields_surface -------------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::register_fields_surface()
{
  // declare and put surface nodal field on part
  nodeSurfaceField_ = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "node_surface_field"));
  stk::mesh::put_field(*nodeSurfaceField_, *superSurfacePart_);
 
  // declare and put surface field on part
  stk::topology::rank_t sideRank = static_cast<stk::topology::rank_t>(metaData_->side_rank());
  surfaceField_ = &(metaData_->declare_field<GenericFieldType>(sideRank, "surface_field"));
  stk::mesh::put_field(*surfaceField_, *superSurfacePart_, 8);
}
  
//--------------------------------------------------------------------------
//-------- initialize_fields -----------------------------------------------
//--------------------------------------------------------------------------
void
SuperElement::initialize_fields()
{
  // just check on whether or not the nodes are all here on the superElementPart_; define gold standard for three element quad4 mesh
  const stk::mesh::EntityId goldElemNodalOrder[27] = {1, 2, 4, 8, 10, 11, 17, 14, 19,
                                                      8, 4, 5, 7, 17, 12, 18, 15, 20,
                                                      7, 5, 3, 6, 18, 13, 9, 16, 21};
  const stk::mesh::EntityId goldElemId[3] = {4,5,6};
  int goldElemIdCount = 0;
  int goldElemNodalOrderCount = 0;
  bool testElemIdPassed = true;
  bool testElemPassed = true;

  // define element selector
  stk::mesh::Selector s_elem = metaData_->locally_owned_part()
    & stk::mesh::Selector(*superElementPart_);

  stk::mesh::BucketVector const& elem_buckets =
    bulkData_->get_buckets(stk::topology::ELEMENT_RANK, s_elem );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get elem
      stk::mesh::Entity elem = b[k];
      const stk::mesh::EntityId elemId = bulkData_->identifier(elem);
  
      if ( elemId != goldElemId[goldElemIdCount] ) {
        testElemIdPassed = false;
        if ( verboseOutput_ )
          NaluEnv::self().naluOutputP0() << " elem id......FAILED " << elemId << " " << goldElemId[goldElemIdCount] << std::endl;
      }
      else {
        if ( verboseOutput_ )
          NaluEnv::self().naluOutputP0() << " elem id......PASSED " << elemId << " " << goldElemId[goldElemIdCount] << std::endl;
      }
      goldElemIdCount++;

      // number of nodes
      int num_nodes = b.num_nodes(k);

      // relations
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      
      if ( verboseOutput_ )
        NaluEnv::self().naluOutputP0() << "... number of nodes: " << num_nodes
                                       << " for element " << bulkData_->identifier(elem) << std::endl;

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // extract nodes
        double * nodalCoords = stk::mesh::field_data(*coordinates_, node);
      
        if ( verboseOutput_ )
          NaluEnv::self().naluOutputP0() << "Node id: " << bulkData_->identifier(node);
        EXPECT_EQ(goldElemNodalOrder[goldElemNodalOrderCount], bulkData_->identifier(node));
        if ( bulkData_->identifier(node) == goldElemNodalOrder[goldElemNodalOrderCount]) {
          if ( verboseOutput_ )
            NaluEnv::self().naluOutputP0() << " ......PASSED" << std::endl;
        }
        else {
          if ( verboseOutput_ )
            NaluEnv::self().naluOutputP0() << " ......FAILED" << std::endl;
          testElemPassed = false;
        }
        
        if ( verboseOutput_ ) {
          for ( int j = 0; j < nDim_; ++j )
            NaluEnv::self().naluOutputP0() << "     coords[" << j << "] " << nodalCoords[j] << std::endl;
        }
        
        // increment count
        goldElemNodalOrderCount++;
      }
    }
  }

  if ( testElemIdPassed )
    NaluEnv::self().naluOutputP0() << "Element Ids Test            PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Element Ids Test            FAILED" << std::endl;

  if ( testElemPassed )
    NaluEnv::self().naluOutputP0() << "Element Connectivities Test PASSED" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Element Connectivities Test FAILED" << std::endl;
    
  // now check nodes in the mesh based on super element part (same selector as above)
  size_t totalNumNodes = 0;
  size_t goldTotalNumNodes = 21;
  stk::mesh::EntityId goldNodalOrder[21] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                            20, 21, 9, 4, 5, 7, 8, 2, 1, 3, 6};
 
  int goldNodalOrderCount = 0;
  bool testNodalPassed = true;
  
  // define node selector
  stk::mesh::Selector s_node = metaData_->locally_owned_part()
    & stk::mesh::Selector(*superElementPart_);

  stk::mesh::BucketVector const& node_buckets =
  bulkData_->get_buckets(stk::topology::NODE_RANK, s_node );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
       ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    totalNumNodes += length;
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // get node
      stk::mesh::Entity node = b[k];
      
      if ( verboseOutput_ )
        NaluEnv::self().naluOutputP0() << "node identifier: " << bulkData_->identifier(node) 
                                       << " " << goldNodalOrder[goldNodalOrderCount] << std::endl;
      
      EXPECT_EQ(goldNodalOrder[goldNodalOrderCount], bulkData_->identifier(node));
      if ( bulkData_->identifier(node) == goldNodalOrder[goldNodalOrderCount] ) {
        if ( verboseOutput_ )
          NaluEnv::self().naluOutputP0() << " ......PASSED" << std::endl;
      }
      else {
        if (verboseOutput_ )
          NaluEnv::self().naluOutputP0() << " ......FAILED" << std::endl;
        testNodalPassed = false;
      }
      
      // increment count
      goldNodalOrderCount++;
    }
  }
  
  if ( totalNumNodes != goldTotalNumNodes )
    testNodalPassed = false;
  
  if ( testNodalPassed )
    NaluEnv::self().naluOutputP0() << "Nodal iteration Test        PASSED" << std::endl;
  else 
    NaluEnv::self().naluOutputP0() << "Nodal iteration Test        FAILED" << std::endl;
  NaluEnv::self().naluOutputP0() << std::endl;
}

} // namespace naluUnit
} // namespace Sierra

TEST(NaluUnit, SuperElement)
{
    sierra::naluUnit::SuperElement superElement;
    superElement.execute();
}

