/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <ScratchViews.h>

namespace sierra {
namespace nalu {

inline
void gather_elem_node_field(const stk::mesh::FieldBase& field,
                            int numNodes,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double*>& shmemView)
{
  for(int i=0; i<numNodes; ++i) {
    shmemView[i] = *static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
  }
}

inline
void gather_elem_node_tensor_field(const stk::mesh::FieldBase& field,
                            int numNodes,
                            int tensorDim1,
                            int tensorDim2,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double***>& shmemView)
{
  for(int i=0; i<numNodes; ++i) {
    const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
    unsigned counter = 0;
    for(int d1=0; d1<tensorDim1; ++d1) { 
      for(int d2=0; d2<tensorDim2; ++d2) {
        shmemView(i,d1,d2) = dataPtr[counter++];
      }
    }
  }
}

inline
void gather_elem_tensor_field(const stk::mesh::FieldBase& field,
                              stk::mesh::Entity elem,
                              int tensorDim1,
                              int tensorDim2,
                              SharedMemView<double**>& shmemView)
{
  const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elem));
  unsigned counter = 0;
  for(int d1=0; d1<tensorDim1; ++d1) { 
    for(int d2=0; d2<tensorDim2; ++d2) {
      shmemView(d1,d2) = dataPtr[counter++];
    }
  }
}

inline
void gather_elem_node_field_3D(const stk::mesh::FieldBase& field,
                            int numNodes,
                            const stk::mesh::Entity* elemNodes,
                            SharedMemView<double**>& shmemView)
{
  for(int i=0; i<numNodes; ++i) {
    const double* dataPtr = static_cast<const double*>(stk::mesh::field_data(field, elemNodes[i]));
    shmemView(i,0) = dataPtr[0];
    shmemView(i,1) = dataPtr[1];
    shmemView(i,2) = dataPtr[2];
  }
}

inline
void gather_elem_node_field(const stk::mesh::FieldBase& field,
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

int
MasterElementViews::create_master_element_views(
  const TeamHandleType& team,
  const std::set<ELEM_DATA_NEEDED>& dataEnums,
  int nDim, int nodesPerElem,
  int numScsIp, int numScvIp)
{
  int numScalars = 0;
  bool needDeriv = false;
  bool needDetj = false;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case SCS_AREAV:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_AREAV is requested.");
         scs_areav = get_shmem_view_2D(team, numScsIp, nDim);
         numScalars += numScsIp * nDim;
         break;

      case SCS_GRAD_OP:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_GRAD_OP is requested.");
         dndx = get_shmem_view_3D(team, numScsIp, nodesPerElem, nDim);
         numScalars += nodesPerElem * numScsIp * nDim;
         needDeriv = true;
         needDetj = true;
         break;

      case SCS_SHIFTED_GRAD_OP:
        ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_SHIFTED_GRAD_OP is requested.");
        dndx_shifted = get_shmem_view_3D(team, numScsIp, nodesPerElem, nDim);
        numScalars += nodesPerElem * numScsIp * nDim;
        needDeriv = true;
        needDetj = true;
        break;

      case SCS_GIJ:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_GIJ is requested.");
         gijUpper = get_shmem_view_3D(team, numScsIp, nDim, nDim);
         gijLower = get_shmem_view_3D(team, numScsIp, nDim, nDim);
         numScalars += 2 * numScsIp * nDim * nDim;
         needDeriv = true;
         break;

      case SCV_VOLUME:
         ThrowRequireMsg(numScvIp > 0, "ERROR, meSCV must be non-null if SCV_VOLUME is requested.");
         scv_volume = get_shmem_view_1D(team, numScvIp);
         numScalars += numScvIp;
         break;

      default: break;
    }
  }
  if (needDeriv) {
    deriv = get_shmem_view_1D(team, numScsIp*nodesPerElem*nDim);
    numScalars += numScsIp * nodesPerElem * nDim;
  }

  if (needDetj) {
    det_j = get_shmem_view_1D(team, numScsIp);
    numScalars += numScsIp;
  }

  return numScalars;
}

void
MasterElementViews::fill_master_element_views(
  const std::set<ELEM_DATA_NEEDED>& dataEnums,
  SharedMemView<double**>* coordsView,
  MasterElement* meSCS,
  MasterElement* meSCV)
{
  double error = 0.0;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case SCS_AREAV:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_AREAV is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_AREAV requested.");
         meSCS->determinant(1, &((*coordsView)(0,0)), &scs_areav(0,0), &error);
         break;
      case SCS_GRAD_OP:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
         meSCS->grad_op(1, &((*coordsView)(0,0)), &dndx(0,0,0), &deriv(0), &det_j(0), &error);
         break;
      case SCS_SHIFTED_GRAD_OP:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        meSCS->shifted_grad_op(1, &((*coordsView)(0,0)), &dndx_shifted(0,0,0), &deriv(0), &det_j(0), &error);
        break;
      case SCS_GIJ:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GIJ is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GIJ requested.");
         meSCS->gij(&((*coordsView)(0,0)), &gijUpper(0,0,0), &gijLower(0,0,0), &deriv(0));
         break;
      case SCV_VOLUME:
         ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_VOLUME is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCV_VOLUME requested.");
         meSCV->determinant(1, &((*coordsView)(0,0)), &scv_volume(0), &error);
         break;
      default: break;
    }
  }
}

ScratchViews::ScratchViews(const TeamHandleType& team,
             const stk::mesh::BulkData& bulkData,
             stk::topology topo,
             ElemDataRequests& dataNeeded)
  :
    meViews(MAX_COORDS_TYPES),
    hasCoordField(MAX_COORDS_TYPES,false)
{
  /* master elements are allowed to be null if they are not required */
  MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();

  int nDim = bulkData.mesh_meta_data().spatial_dimension();
  int nodesPerElem = topo.num_nodes();
  int numScsIp = meSCS != nullptr ? meSCS->numIntPoints_ : 0;
  int numScvIp = meSCV != nullptr ? meSCV->numIntPoints_ : 0;

  create_needed_field_views(team, dataNeeded, bulkData, nodesPerElem);

  create_needed_master_element_views(team, dataNeeded, nDim, nodesPerElem, numScsIp, numScvIp);
}

void ScratchViews::create_needed_field_views(const TeamHandleType& team,
                               const ElemDataRequests& dataNeeded,
                               const stk::mesh::BulkData& bulkData,
                               int nodesPerElem)
{
  int numScalars = 0;
  const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
  unsigned numFields = meta.get_fields().size();
  fieldViews.resize(numFields, nullptr);

  const FieldSet& neededFields = dataNeeded.get_fields();
  for(const FieldInfo& fieldInfo : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
    ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
    unsigned scalarsDim1 = fieldInfo.scalarsDim1;
    unsigned scalarsDim2 = fieldInfo.scalarsDim2;

    if (fieldEntityRank==stk::topology::ELEM_RANK) {
      if (scalarsDim2 == 0) {
        fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double*>>(get_shmem_view_1D(team, scalarsDim1));
        numScalars += scalarsDim1;
      }
      else {
        fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double**>>(get_shmem_view_2D(team, scalarsDim1, scalarsDim2));
        numScalars += scalarsDim1 * scalarsDim2;
      }
    }
    else if (fieldEntityRank==stk::topology::NODE_RANK) {
      if (scalarsDim2 == 0) {
        if (scalarsDim1 == 1) {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double*>>(get_shmem_view_1D(team, nodesPerElem));
          numScalars += nodesPerElem;
        }
        else {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double**>>(get_shmem_view_2D(team, nodesPerElem, scalarsDim1));
          numScalars += nodesPerElem * scalarsDim1;
        }
      }
      else {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<double***>>(get_shmem_view_3D(team, nodesPerElem, scalarsDim1, scalarsDim2));
          numScalars += nodesPerElem * scalarsDim1 * scalarsDim2;
      }
    }
    else {
      ThrowRequireMsg(false,"Only elem-rank and node-rank fields supported for scratch-views currently.");
    }
  }

  // Track total bytes required for field allocations
  num_bytes_required += numScalars * sizeof(double);
}

void ScratchViews::create_needed_master_element_views(const TeamHandleType& team,
                                        const ElemDataRequests& dataNeeded,
                                        int nDim, int nodesPerElem,
                                        int numScsIp, int numScvIp)
{
  int numScalars = 0;

  for (auto it = dataNeeded.get_coordinates_map().begin();
       it != dataNeeded.get_coordinates_map().end(); ++it) {
    hasCoordField[it->first] = true;
    numScalars += meViews[it->first].create_master_element_views(
      team, dataNeeded.get_data_enums(it->first),
      nDim, nodesPerElem, numScsIp, numScvIp);
  }

  num_bytes_required += numScalars * sizeof(double);
}

int get_num_bytes_pre_req_data(ElemDataRequests& dataNeededBySuppAlgs, int nDim)
{
  /* master elements are allowed to be null if they are not required */
  MasterElement *meSCS = dataNeededBySuppAlgs.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeededBySuppAlgs.get_cvfem_volume_me();
  
  const int nodesPerElem = meSCS != nullptr ? meSCS->nodesPerElement_ :
    (meSCV != nullptr ? meSCV->nodesPerElement_ : 0);
  const int numScsIp = meSCS != nullptr ? meSCS->numIntPoints_ : 0;
  const int numScvIp = meSCV != nullptr ? meSCV->numIntPoints_ : 0;
  int numBytes = 0;
  
  const FieldSet& neededFields = dataNeededBySuppAlgs.get_fields();
  for(const FieldInfo& fieldInfo : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
    ThrowAssertMsg(fieldEntityRank == stk::topology::NODE_RANK || fieldEntityRank == stk::topology::ELEM_RANK, "Currently only node and element fields are supported.");
    unsigned scalarsPerEntity = fieldInfo.scalarsDim1;
    unsigned entitiesPerElem = fieldEntityRank==stk::topology::ELEM_RANK ? 1 : nodesPerElem;

    // Catch errors if user requests nodal field but has not registered any
    // MasterElement we need to get nodesPerElem
    ThrowRequire(entitiesPerElem > 0);
    if (fieldInfo.scalarsDim2 > 1) {
      scalarsPerEntity *= fieldInfo.scalarsDim2;
    }
    numBytes += entitiesPerElem*scalarsPerEntity*sizeof(double);
  }
  
  for (auto it = dataNeededBySuppAlgs.get_coordinates_map().begin();
       it != dataNeededBySuppAlgs.get_coordinates_map().end(); ++it)
  {
    const std::set<ELEM_DATA_NEEDED>& dataEnums =
      dataNeededBySuppAlgs.get_data_enums(it->first);
    int dndxLength = 0, gUpperLength = 0, gLowerLength = 0;

    // Updated logic for data sharing of deriv and det_j
    bool needDeriv = false;
    bool needDetj = false;
    for(ELEM_DATA_NEEDED data : dataEnums) {
      switch(data)
      {
      case SCS_AREAV: numBytes += nDim * numScsIp * sizeof(double);
        break;
      case SCS_GRAD_OP:
      case SCS_SHIFTED_GRAD_OP:
        dndxLength = nodesPerElem*numScsIp*nDim;
        needDeriv = true;
        needDetj = true;
        numBytes += dndxLength * sizeof(double);
        break;
      case SCV_VOLUME: numBytes += numScvIp * sizeof(double);
        break;
      case SCS_GIJ: 
        gUpperLength = nDim*nDim*numScsIp;
        gLowerLength = nDim*nDim*numScsIp;
        needDeriv = true;
        numBytes += (gUpperLength + gLowerLength ) * sizeof(double);
        break;
      default: break;
      }
    }

    if (needDeriv)
      numBytes += nodesPerElem*numScsIp*nDim * sizeof(double);

    if (needDetj)
      numBytes += numScsIp * sizeof(double);
  }

  // Add a 64 byte padding to the buffer size requested
  return numBytes + 8 * sizeof(double);
}

void fill_pre_req_data(
  ElemDataRequests& dataNeeded,
  const stk::mesh::BulkData& bulkData,
  stk::topology topo,
  stk::mesh::Entity elem,
  ScratchViews& prereqData)
{
  int nodesPerElem = topo.num_nodes();

  MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();
  prereqData.elemNodes = bulkData.begin_nodes(elem);

  const FieldSet& neededFields = dataNeeded.get_fields();
  for(const FieldInfo& fieldInfo : neededFields) {
    stk::mesh::EntityRank fieldEntityRank = fieldInfo.field->entity_rank();
    unsigned scalarsDim1 = fieldInfo.scalarsDim1;
    bool isTensorField = fieldInfo.scalarsDim2 > 1;

    if (fieldEntityRank==stk::topology::ELEM_RANK) {
      if (isTensorField) {
        SharedMemView<double**>& shmemView = prereqData.get_scratch_view_2D(*fieldInfo.field);
        gather_elem_tensor_field(*fieldInfo.field, elem, scalarsDim1, fieldInfo.scalarsDim2, shmemView);
      }
      else {
        SharedMemView<double*>& shmemView = prereqData.get_scratch_view_1D(*fieldInfo.field);
        unsigned len = shmemView.dimension(0);
        double* fieldDataPtr = static_cast<double*>(stk::mesh::field_data(*fieldInfo.field, elem));
        for(unsigned i=0; i<len; ++i) {
          shmemView(i) = fieldDataPtr[i];
        }
      }
    }
    else if (fieldEntityRank == stk::topology::NODE_RANK) {
      if (isTensorField) {
        SharedMemView<double***>& shmemView3D = prereqData.get_scratch_view_3D(*fieldInfo.field);
        gather_elem_node_tensor_field(*fieldInfo.field, nodesPerElem, scalarsDim1, fieldInfo.scalarsDim2, bulkData.begin_nodes(elem), shmemView3D);
      }
      else {
        if (scalarsDim1 == 1) {
          SharedMemView<double*>& shmemView1D = prereqData.get_scratch_view_1D(*fieldInfo.field);
          gather_elem_node_field(*fieldInfo.field, nodesPerElem, prereqData.elemNodes, shmemView1D);
        }
        else {
          SharedMemView<double**>& shmemView2D = prereqData.get_scratch_view_2D(*fieldInfo.field);
          if (scalarsDim1 == 3) {
            gather_elem_node_field_3D(*fieldInfo.field, nodesPerElem, prereqData.elemNodes, shmemView2D);
          }
          else {
            gather_elem_node_field(*fieldInfo.field, nodesPerElem, scalarsDim1, prereqData.elemNodes, shmemView2D);
          }
        }
      }
    }
    else {
      ThrowRequireMsg(false, "Only node and element fields supported currently.");
    }
  } 


  for (auto it = dataNeeded.get_coordinates_map().begin();
       it != dataNeeded.get_coordinates_map().end(); ++it) {
    auto cType = it->first;
    auto coordField = it->second;

    const std::set<ELEM_DATA_NEEDED>& dataEnums = dataNeeded.get_data_enums(cType);
    SharedMemView<double**>* coordsView = &prereqData.get_scratch_view_2D(*coordField);
    MasterElementViews& meData = prereqData.get_me_views(cType);

    meData.fill_master_element_views(dataEnums, coordsView, meSCS, meSCV);
  }
}

}
}

