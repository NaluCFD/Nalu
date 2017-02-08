/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScratchViews_h
#define ScratchViews_h

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <ElemDataRequests.h>

#include <KokkosInterface.h>

#include <set>

namespace sierra{
namespace nalu{

inline
void gather_elem_node_field(const stk::mesh::FieldBase& field,
                            stk::mesh::Entity elem,
                            int numNodes,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double*>& shmemView)
{
  for(int i=0; i<numNodes; ++i) {
    const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
    shmemView[i] = *dataPtr;
  }
}

inline
void gather_elem_node_field(const stk::mesh::FieldBase& field,
                            stk::mesh::Entity elem,
                            int numNodes,
                            int scalarsPerNode,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double**>& shmemView)
{
  for(int i=0; i<numNodes; ++i) {
    const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
    for(int d=0; d<scalarsPerNode; ++d) {
      shmemView(i,d) = dataPtr[d];
    }    
  }
}

class ScratchViews
{
public:
  ScratchViews(const TeamHandleType& team,
               const stk::mesh::BulkData& bulkData,
               stk::topology topo,
               ElemDataRequests& dataNeeded)
  : elemNodes(), scs_areav(), dndx(), deriv(), det_j()
  {
    // master elements are allowed to be null if they are not required
    MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
    MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();

    int nDim = bulkData.mesh_meta_data().spatial_dimension();
    int nodesPerElem = topo.num_nodes();
    int numScsIp = meSCS != nullptr ? meSCS->numIntPoints_ : 0;
    int numScvIp = meSCV != nullptr ? meSCV->numIntPoints_ : 0;

    create_needed_field_views(team, dataNeeded, bulkData, nodesPerElem);

    create_needed_master_element_views(team, dataNeeded, nDim, nodesPerElem, numScsIp, numScvIp);
  }

  SharedMemView<double*>& get_scratch_view_1D(const stk::mesh::FieldBase& field)
  {
    return fieldViews1D[field.mesh_meta_data_ordinal()];
  }

  SharedMemView<double**>& get_scratch_view_2D(const stk::mesh::FieldBase& field)
  {
    return fieldViews2D[field.mesh_meta_data_ordinal()];
  }

  SharedMemView<double***>& get_scratch_view_3D(const stk::mesh::FieldBase& field)
  {
    return fieldViews3D[field.mesh_meta_data_ordinal()];
  }

  SharedMemView<stk::mesh::Entity*> elemNodes;
  SharedMemView<double**> scs_areav;
  SharedMemView<double***> dndx;
  SharedMemView<double*> deriv;
  SharedMemView<double*> det_j;
  SharedMemView<double*> scv_volume;

private:
  void create_needed_field_views(const TeamHandleType& team,
                                 const ElemDataRequests& dataNeeded,
                                 const stk::mesh::BulkData& bulkData,
                                 int nodesPerElem)
  {
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    unsigned numFields = meta.get_fields().size();
    // FIXME: Ideally, size based on what was actually added; waht about fast access?
    fieldViews1D.resize(numFields);
    fieldViews2D.resize(numFields);
    fieldViews3D.resize(numFields);

    const FieldSet& neededFields = dataNeeded.get_fields();
    for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
      stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
      ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
      unsigned scalarsPerEntity = fieldPtr->max_size(fieldEntityRank);
      if (fieldEntityRank==stk::topology::ELEM_RANK) {
        fieldViews1D[fieldPtr->mesh_meta_data_ordinal()] = get_shmem_view_1D(team, scalarsPerEntity);
      }
      else if (fieldEntityRank==stk::topology::NODE_RANK) {
        if (scalarsPerEntity == 1) {
          fieldViews1D[fieldPtr->mesh_meta_data_ordinal()] = get_shmem_view_1D(team, nodesPerElem);
        }
        else {
          fieldViews2D[fieldPtr->mesh_meta_data_ordinal()] = get_shmem_view_2D(team, nodesPerElem, scalarsPerEntity);
        }
      }
      else {
        ThrowRequireMsg(false,"Only elem-rank and node-rank fields supported for scratch-views currently.");
      }
    }
  }

  void create_needed_master_element_views(const TeamHandleType& team,
                                          const ElemDataRequests& dataNeeded,
                                          int nDim, int nodesPerElem,
                                          int numScsIp, int numScvIp)
  {
    const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums();
    for(ELEM_DATA_NEEDED data : dataEnums) {
      switch(data)
      {
        case NODES:
           elemNodes = get_entity_shmem_view_1D(team, nodesPerElem);
           break;

        case SCS_AREAV:
           ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_AREAV is requested.");
           scs_areav = get_shmem_view_2D(team, numScsIp, nDim);
           break;

        case SCS_GRAD_OP:
           ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_GRAD_OP is requested.");
           dndx = get_shmem_view_3D(team, numScsIp, nodesPerElem, nDim);
           deriv = get_shmem_view_1D(team, numScsIp*nodesPerElem*nDim);
           det_j = get_shmem_view_1D(team, numScsIp);
           break;

        case SCV_VOLUME:
           ThrowRequireMsg(numScvIp > 0, "ERROR, meSCV must be non-null if SCV_VOLUME is requested.");
           scv_volume = get_shmem_view_1D(team, numScvIp);
           break;

        default: break;
      }
    }
  }

  std::vector<SharedMemView<double*> > fieldViews1D;
  std::vector<SharedMemView<double**> > fieldViews2D;
  std::vector<SharedMemView<double***> > fieldViews3D;
};

inline
int get_num_bytes_pre_req_data(
  ElemDataRequests& dataNeededBySuppAlgs,
  int nDim)
{
  // master elements are allowed to be null if they are not required
  MasterElement *meSCS = dataNeededBySuppAlgs.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeededBySuppAlgs.get_cvfem_volume_me();
  
  const int nodesPerElem = meSCS != nullptr ? meSCS->nodesPerElement_ : 0;
  const int numScsIp = meSCS != nullptr ? meSCS->numIntPoints_ : 0;
  const int numScvIp = meSCV != nullptr ? meSCV->numIntPoints_ : 0;
  int numBytes = 0;
  
  const FieldSet& neededFields = dataNeededBySuppAlgs.get_fields();
  for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
    ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
    unsigned scalarsPerEntity = fieldPtr->max_size(fieldEntityRank);
    unsigned entitiesPerElem = fieldEntityRank==stk::topology::ELEM_RANK ? 1 : nodesPerElem;
    numBytes += entitiesPerElem*scalarsPerEntity*sizeof(double);
  }
  
  const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeededBySuppAlgs.get_data_enums();
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
      {
      case NODES: numBytes += nodesPerElem * sizeof(stk::mesh::Entity);
        break;
      case SCS_AREAV: numBytes += nDim * numScsIp * sizeof(double);
        break;
      case SCS_GRAD_OP: numBytes += (nodesPerElem*numScsIp*nDim*2 + numScsIp) * sizeof(double); // dndx, deriv and detJ
        break;
      case SCV_VOLUME: numBytes += numScvIp * sizeof(double);
      default: break;
      }
  }
  
  return numBytes*2;
}

inline
void fill_pre_req_data(
  ElemDataRequests& dataNeeded,
  const stk::mesh::BulkData& bulkData,
  stk::topology topo,
  stk::mesh::Entity elem,
  const stk::mesh::FieldBase* coordField,
  ScratchViews& prereqData)
{
  int nodesPerElem = topo.num_nodes();
  const stk::mesh::Entity* nodes = nullptr;

  // extract master elements
  MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();

  const FieldSet& neededFields = dataNeeded.get_fields();
  for(const stk::mesh::FieldBase* fieldPtr : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldPtr->entity_rank();
    if (fieldEntityRank==stk::topology::ELEM_RANK) {
      SharedMemView<double*>& shmemView = prereqData.get_scratch_view_1D(*fieldPtr);
      unsigned len = shmemView.dimension(0);
      double* fieldDataPtr = static_cast<double*>(stk::mesh::field_data(*fieldPtr, elem));
      for(unsigned i=0; i<len; ++i) {
        shmemView(i) = fieldDataPtr[i];
      }
    }
    else if (fieldEntityRank == stk::topology::NODE_RANK) {
      SharedMemView<double*>& shmemView1D = prereqData.get_scratch_view_1D(*fieldPtr);
      if (shmemView1D.dimension(0) > 0) {
        gather_elem_node_field(*fieldPtr, elem, nodesPerElem, bulkData.begin_nodes(elem), shmemView1D);
      }
      else {
        SharedMemView<double**>& shmemView2D = prereqData.get_scratch_view_2D(*fieldPtr);
        int scalarsPerNode = shmemView2D.dimension(1);
        gather_elem_node_field(*fieldPtr, elem, nodesPerElem, scalarsPerNode, bulkData.begin_nodes(elem), shmemView2D);
      }
    }
    else {
      ThrowRequireMsg(false, "Only node and element fields supported currently.");
    }
  }

  SharedMemView<double**>& coords = prereqData.get_scratch_view_2D(*coordField);
  const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums();
  double error = 0;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case NODES:
         nodes = bulkData.begin_nodes(elem);
         for(int i=0; i<nodesPerElem; ++i) {
           prereqData.elemNodes(i) = nodes[i];
         }
         break;

      case SCS_AREAV:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_AREAV is requested.");
         meSCS->determinant(1, &coords(0,0), &prereqData.scs_areav(0,0), &error);
         break;

      case SCS_GRAD_OP:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
         meSCS->grad_op(1, &coords(0,0), &prereqData.dndx(0,0,0), &prereqData.deriv(0), &prereqData.det_j(0), &error);
         break;
      case SCV_VOLUME:
         ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_VOLUME is requested.");
         meSCV->determinant(1, &coords(0,0), &prereqData.scv_volume(0), &error);

      default: break;
    }
  }
}

} // namespace nalu
} // namespace Sierra

#endif
