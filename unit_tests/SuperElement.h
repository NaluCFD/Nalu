/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SuperElement_h
#define SuperElement_h

// stk_mesh
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

// stk_search
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>

// STL
#include <vector>
#include <map>

// field types
typedef stk::mesh::Field<double>  ScalarFieldType;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorFieldType;
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;

namespace stk {
  namespace io {
    class StkMeshIoBroker;
  }
  namespace mesh {
    class Part;
    class MetaData;
    class BulkData;
    typedef std::vector<Part*> PartVector;
    typedef std::vector<stk::mesh::EntityId> EntityIdVector;
    struct Entity;
  }
}

namespace sierra {
namespace naluUnit {

struct EntityNodeSharing {
public:
  stk::mesh::EntityId edgeId_;
  stk::mesh::EntityId globalNodeId_;
  int localIndex_;
};

class SuperElement
{
public:

  // constructor/destructor
  SuperElement();
  ~SuperElement();

  void execute();

  // part for super part
  void declare_super_part();
  void declare_super_part_surface();

  // part for edges and faces
  void declare_edge_part();
  void declare_face_part();

  // find size of each entity on the mesh 
  void size_of_edges();
  void size_of_faces();
  void size_of_elements();

  void create_nodes();

  void consolidate_node_ids();
  
  void create_elements();
  void create_elements_surface();

  // create and destroy faces off of original part
  void create_edges();
  void create_faces();
  void delete_edges();
  void delete_faces();
  
  // register nodal and elemental fields
  void register_fields();
  void register_fields_surface();
  
  // initialize nodal fields
  void initialize_fields();

  void initialize_mesh(stk::ParallelMachine pm);
  void createThreeElemQuadMesh_1proc(stk::mesh::BulkData& mesh);

  const int pOrder_;
  
  // aura on/off
  const bool activateAura_;

  double currentTime_;
  size_t resultsFileIndex_;
  int nDim_;

  // meta, bulk and io
  stk::mesh::MetaData *metaData_;
  stk::mesh::BulkData *bulkData_;
  stk::io::StkMeshIoBroker *ioBroker_;

  // fields
  ScalarFieldType *nodeField_;
  GenericFieldType *elemField_;
  ScalarFieldType *nodeSurfaceField_;
  GenericFieldType *surfaceField_;
  VectorFieldType *coordinates_;

  // the original [volume] part of the lower order mesh, e.g., block_1, Topo::quad4
  std::string originalPartName_; 
  std::string originalSurfacePartName_;
  
  // the new volume part for the higher order mesh, e.g., block_1_SE, Topo::superElement
  std::string superElementPartName_;
  std::string superElementSurfacePartName_;

  // the set of nodes that are promoted
  std::string promotedNodesPartName_;
  
  // name for the new edges and faces created
  std::string edgePartName_;
  std::string facePartName_;

  // verbosity level for output
  const bool verboseOutput_;

  // part associated with lower order standard element
  stk::mesh::Part *originalPart_;
  stk::mesh::Part *originalSurfacePart_;

  // part associated with super element
  stk::mesh::Part *superElementPart_;
  stk::mesh::Part *superSurfacePart_;

  // in-transit part associated with augmented/promoted nodes
  stk::mesh::Part *promotedNodesPart_;

  // in-transit part associated with edges and faces
  stk::mesh::Part *edgePart_;
  stk::mesh::Part *facePart_;

  // entity count to define node count
  size_t numberOfEdges_;
  size_t numberOfFaces_;
  size_t numberOfElements_;

  // vector of new nodes
  std::vector<stk::mesh::Entity> promotedNodesVec_;

  // create mapping of parent ids nodes to the new node
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> > parentEdgeNodesMap_;
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> > parentFaceNodesMap_;
  std::map<stk::mesh::EntityId, std::vector<stk::mesh::Entity> > parentElemNodesMap_;
  
  // mapping of low order element id to the new higher order element
  std::map<stk::mesh::EntityId, stk::mesh::Entity > superElementElemMap_;
  
  // keep something that holds the set of edges and who besides local rank holds them
  std::vector<std::vector<int> > sharedProcsEdge_;

  // communication patter for edge-nodes
  std::map<int, std::vector<EntityNodeSharing> > edgeNodeSharingMap_;
  std::map<int, std::vector<EntityNodeSharing> > edgeNodeSharingOffProcMap_;
};

} // namespace naluUnit
} // namespace Sierra

#endif
