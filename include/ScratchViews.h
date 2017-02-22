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

struct ViewHolder {
  virtual ~ViewHolder() {}
};

template<typename T>
struct ViewT : public ViewHolder {
  ViewT(T view) : view_(view) {}
  virtual ~ViewT(){}
  T view_;
};

class ScratchViews
{
public:
  ScratchViews(const TeamHandleType& team,
               const stk::mesh::BulkData& bulkData,
               stk::topology topo,
               ElemDataRequests& dataNeeded)
  : elemNodes(nullptr), scs_areav(), dndx(), deriv(), det_j()
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

  virtual ~ScratchViews() {
    for(ViewHolder* vh : fieldViews) {
      delete vh;
    }
  }

  SharedMemView<double*>& get_scratch_view_1D(const stk::mesh::FieldBase& field)
  {
    ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 1D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
    ViewT<SharedMemView<double*>>* vt = static_cast<ViewT<SharedMemView<double*>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
    return vt->view_;
  }

  SharedMemView<double**>& get_scratch_view_2D(const stk::mesh::FieldBase& field)
  {
    ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 2D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
    ViewT<SharedMemView<double**>>* vt = static_cast<ViewT<SharedMemView<double**>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
    return vt->view_;
  }

  SharedMemView<double***>& get_scratch_view_3D(const stk::mesh::FieldBase& field)
  {
    ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 3D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
    ViewT<SharedMemView<double***>>* vt = static_cast<ViewT<SharedMemView<double***>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
    return vt->view_;
  }

  SharedMemView<double****>& get_scratch_view_4D(const stk::mesh::FieldBase& field)
  {
    ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 4D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
    ViewT<SharedMemView<double****>>* vt = static_cast<ViewT<SharedMemView<double****>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
    return vt->view_;
  }

  const stk::mesh::Entity* elemNodes;
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
    fieldViews.resize(numFields, nullptr);

    const FieldSet& neededFields = dataNeeded.get_fields();
    for(const FieldInfo& fieldInfo : neededFields) {
      stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
      ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
      unsigned scalarsPerEntity = fieldInfo.scalarsPerEntity;
      if (fieldEntityRank==stk::topology::ELEM_RANK) {
        fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double*>>(get_shmem_view_1D(team, scalarsPerEntity));
      }
      else if (fieldEntityRank==stk::topology::NODE_RANK) {
        if (scalarsPerEntity == 1) {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double*>>(get_shmem_view_1D(team, nodesPerElem));
        }
        else {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double**>>(get_shmem_view_2D(team, nodesPerElem, scalarsPerEntity));
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

  std::vector<ViewHolder*> fieldViews;
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
  for(const FieldInfo& fieldInfo : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
    ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
    unsigned scalarsPerEntity = fieldInfo.scalarsPerEntity;
    unsigned entitiesPerElem = fieldEntityRank==stk::topology::ELEM_RANK ? 1 : nodesPerElem;
    numBytes += entitiesPerElem*scalarsPerEntity*sizeof(double);
  }
  
  const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeededBySuppAlgs.get_data_enums();
  int dndxLength = 0, derivLength = 0, detJLength = 0;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
      {
      case SCS_AREAV: numBytes += nDim * numScsIp * sizeof(double);
        break;
      case SCS_GRAD_OP:
        dndxLength = nodesPerElem*numScsIp*nDim;
        derivLength = nodesPerElem*numScsIp*nDim;
        detJLength = numScsIp;
        numBytes += (dndxLength + derivLength + detJLength) * sizeof(double);
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

  // extract master elements
  MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();

  const FieldSet& neededFields = dataNeeded.get_fields();
  for(const FieldInfo& fieldInfo : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
    unsigned scalarsPerEntity = fieldInfo.scalarsPerEntity;
    if (fieldEntityRank==stk::topology::ELEM_RANK) {
      SharedMemView<double*>& shmemView = prereqData.get_scratch_view_1D(*fieldInfo.field);
      unsigned len = shmemView.dimension(0);
      double* fieldDataPtr = static_cast<double*>(stk::mesh::field_data(*fieldInfo.field, elem));
      for(unsigned i=0; i<len; ++i) {
        shmemView(i) = fieldDataPtr[i];
      }
    }
    else if (fieldEntityRank == stk::topology::NODE_RANK) {
      if (scalarsPerEntity == 1) {
        SharedMemView<double*>& shmemView1D = prereqData.get_scratch_view_1D(*fieldInfo.field);
        gather_elem_node_field(*fieldInfo.field, elem, nodesPerElem, bulkData.begin_nodes(elem), shmemView1D);
      }
      else {
        SharedMemView<double**>& shmemView2D = prereqData.get_scratch_view_2D(*fieldInfo.field);
        gather_elem_node_field(*fieldInfo.field, elem, nodesPerElem, scalarsPerEntity, bulkData.begin_nodes(elem), shmemView2D);
      }
    }
    else {
      ThrowRequireMsg(false, "Only node and element fields supported currently.");
    }
  }

  SharedMemView<double**>& coords = prereqData.get_scratch_view_2D(*coordField);
  prereqData.elemNodes = bulkData.begin_nodes(elem);
  const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums();
  double error = 0;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
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
