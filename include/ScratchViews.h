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
#include <master_element/MasterElement.h>
#include <KokkosInterface.h>

#include <set>
#include <type_traits>

namespace sierra{
namespace nalu{

struct ViewHolder {
  virtual ~ViewHolder() {}
  int dim_;
};

template<typename T>
struct ViewT : public ViewHolder {
  ViewT(T view, int dim) : view_(view) {dim_ = dim;}
  virtual ~ViewT(){}
  T view_;
};

template<typename T>
class MasterElementViews
{
public:
  typedef T value_type;

  MasterElementViews() = default;
  virtual ~MasterElementViews() = default;

  int create_master_element_views(
    const TeamHandleType& team,
    const std::set<ELEM_DATA_NEEDED>& dataEnums,
    int nDim, int nodesPerFace, int nodesPerElem,
    int numFaceIp, int numScsIp, int numScvIp, int numFemIp);

  void fill_master_element_views(
    const std::set<ELEM_DATA_NEEDED>& dataEnums,
    SharedMemView<double**>* coordsView,
    MasterElement* meFC,
    MasterElement* meSCS,
    MasterElement* meSCV,
    MasterElement* meFEM,
    const int* faceOrdinals = nullptr);

  void fill_master_element_views_new_me(
    const std::set<ELEM_DATA_NEEDED>& dataEnums,
    SharedMemView<DoubleType**>* coordsView,
    MasterElement* meFC,
    MasterElement* meSCS,
    MasterElement* meSCV,
    MasterElement* meFEM,
    const int* faceOrdinals = nullptr);

  SharedMemView<T**> fc_areav;
  SharedMemView<T**> scs_areav;
  SharedMemView<T***> dndx_fc_scs;
  SharedMemView<T***> dndx;
  SharedMemView<T***> dndx_shifted;
  SharedMemView<T***> dndx_scv;
  SharedMemView<T***> dndx_fem;
  SharedMemView<T***> deriv_fc_scs;
  SharedMemView<T***> deriv;
  SharedMemView<T***> deriv_scv;
  SharedMemView<T***> deriv_fem;
  SharedMemView<T*> det_j_fc_scs;
  SharedMemView<T*> det_j;
  SharedMemView<T*> det_j_scv;
  SharedMemView<T*> det_j_fem;
  SharedMemView<T*> scv_volume;
  SharedMemView<T***> gijUpper;
  SharedMemView<T***> gijLower;
};

template<typename T>
class ScratchViews
{
public:
  typedef T value_type;

  ScratchViews(const TeamHandleType& team,
               const stk::mesh::BulkData& bulkData,
               int nodesPerEntity,
               const ElemDataRequests& dataNeeded);

  virtual ~ScratchViews() {
    for(ViewHolder* vh : fieldViews) {
      delete vh;
    }
  }

  inline
  SharedMemView<T*>& get_scratch_view_1D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<T**>& get_scratch_view_2D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<T***>& get_scratch_view_3D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<T****>& get_scratch_view_4D(const stk::mesh::FieldBase& field);

  inline
  MasterElementViews<T>& get_me_views(const COORDS_TYPES cType)
  {
    ThrowRequire(hasCoordField[cType] == true);
    return meViews[cType];
  }
  inline bool has_coord_field(const COORDS_TYPES cType) const { return hasCoordField[cType]; }

  inline int total_bytes() const { return num_bytes_required; }

  const stk::mesh::Entity* elemNodes;

  inline const std::vector<ViewHolder*>& get_field_views() const { return fieldViews; }

private:
  void create_needed_field_views(const TeamHandleType& team,
                                 const ElemDataRequests& dataNeeded,
                                 const stk::mesh::BulkData& bulkData,
                                 int nodesPerElem);

  void create_needed_master_element_views(const TeamHandleType& team,
                                          const ElemDataRequests& dataNeeded,
                                          int nDim, int nodesPerFace, int nodesPerElem,
                                          int numFaceIp, int numScsIp, int numScvIp, int numFemIp);

  std::vector<ViewHolder*> fieldViews;
  MasterElementViews<T> meViews[MAX_COORDS_TYPES];
  bool hasCoordField[MAX_COORDS_TYPES] = {false, false};
  int num_bytes_required{0};
};

template<typename T>
SharedMemView<T*>& ScratchViews<T>::get_scratch_view_1D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 1D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<T*>>* vt = static_cast<ViewT<SharedMemView<T*>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

template<typename T>
SharedMemView<T**>& ScratchViews<T>::get_scratch_view_2D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 2D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<T**>>* vt = static_cast<ViewT<SharedMemView<T**>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

template<typename T>
SharedMemView<T***>& ScratchViews<T>::get_scratch_view_3D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 3D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<T***>>* vt = static_cast<ViewT<SharedMemView<T***>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

template<typename T>
SharedMemView<T****>& ScratchViews<T>::get_scratch_view_4D(const stk::mesh::FieldBase& field)
{
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 4D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<T****>>* vt = static_cast<ViewT<SharedMemView<T****>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

template<typename T>
int MasterElementViews<T>::create_master_element_views(
  const TeamHandleType& team,
  const std::set<ELEM_DATA_NEEDED>& dataEnums,
  int nDim, int nodesPerFace, int nodesPerElem,
  int numFaceIp, int numScsIp, int numScvIp, int numFemIp)
{
  int numScalars = 0;
  bool needDeriv = false; bool needDerivScv = false; bool needDerivFem = false; bool needDerivFC = false;
  bool needDetj = false; bool needDetjScv = false; bool needDetjFem = false; bool needDetjFC = false;
  bool femGradOp = false; bool femShiftedGradOp = false;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case FC_AREAV:
          ThrowRequireMsg(numFaceIp > 0, "ERROR, meFC must be non-null if FC_AREAV is requested.");
          fc_areav = get_shmem_view_2D<T>(team, numFaceIp, nDim);
          numScalars += numFaceIp * nDim;
          break;
      case SCS_FACE_GRAD_OP:
          ThrowRequireMsg(numFaceIp > 0, "ERROR, meSCS must be non-null if SCS_FACE_GRAD_OP is requested.");
          dndx_fc_scs = get_shmem_view_3D<T>(team, numFaceIp, nodesPerElem, nDim);
          numScalars += nodesPerElem * numFaceIp * nDim;
          needDerivFC = true;
          needDetjFC = true;
          break;
      case SCS_AREAV:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_AREAV is requested.");
         scs_areav = get_shmem_view_2D<T>(team, numScsIp, nDim);
         numScalars += numScsIp * nDim;
         break;

      case SCS_GRAD_OP:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_GRAD_OP is requested.");
         dndx = get_shmem_view_3D<T>(team, numScsIp, nodesPerElem, nDim);
         numScalars += nodesPerElem * numScsIp * nDim;
         needDeriv = true;
         needDetj = true;
         break;

      case SCS_SHIFTED_GRAD_OP:
        ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_SHIFTED_GRAD_OP is requested.");
        dndx_shifted = get_shmem_view_3D<T>(team, numScsIp, nodesPerElem, nDim);
        numScalars += nodesPerElem * numScsIp * nDim;
        needDeriv = true;
        needDetj = true;
        break;

      case SCS_GIJ:
         ThrowRequireMsg(numScsIp > 0, "ERROR, meSCS must be non-null if SCS_GIJ is requested.");
         gijUpper = get_shmem_view_3D<T>(team, numScsIp, nDim, nDim);
         gijLower = get_shmem_view_3D<T>(team, numScsIp, nDim, nDim);
         numScalars += 2 * numScsIp * nDim * nDim;
         needDeriv = true;
         break;

      case SCV_VOLUME:
         ThrowRequireMsg(numScvIp > 0, "ERROR, meSCV must be non-null if SCV_VOLUME is requested.");
         scv_volume = get_shmem_view_1D<T>(team, numScvIp);
         numScalars += numScvIp;
         break;

      case SCV_GRAD_OP:
         ThrowRequireMsg(numScvIp > 0, "ERROR, meSCV must be non-null if SCV_GRAD_OP is requested.");
         dndx_scv = get_shmem_view_3D<T>(team, numScvIp, nodesPerElem, nDim);
         numScalars += nodesPerElem * numScvIp * nDim;
         needDerivScv = true;
         needDetjScv = true;
         break;

      case FEM_GRAD_OP:
         ThrowRequireMsg(numFemIp > 0, "ERROR, meFEM must be non-null if FEM_GRAD_OP is requested.");
         dndx_fem = get_shmem_view_3D<T>(team, numFemIp, nodesPerElem, nDim);
         numScalars += nodesPerElem * numFemIp * nDim;
         needDerivFem = true;
         needDetjFem = true;
         femGradOp = true;
         break;

      case FEM_SHIFTED_GRAD_OP:
         ThrowRequireMsg(numFemIp > 0, "ERROR, meFEM must be non-null if FEM_SHIFTED_GRAD_OP is requested.");
         dndx_fem = get_shmem_view_3D<T>(team, numFemIp, nodesPerElem, nDim);
         numScalars += nodesPerElem * numFemIp * nDim;
         needDerivFem = true;
         needDetjFem = true;
         femShiftedGradOp = true;
         break;

      default: break;
    }
  }

  if (needDerivFC) {
    deriv_fc_scs = get_shmem_view_3D<T>(team, numFaceIp,nodesPerElem,nDim);
    numScalars += numFaceIp * nodesPerElem * nDim;
  }

  if (needDeriv) {
    deriv = get_shmem_view_3D<T>(team, numScsIp,nodesPerElem,nDim);
    numScalars += numScsIp * nodesPerElem * nDim;
  }

  if (needDerivScv) {
    deriv_scv = get_shmem_view_3D<T>(team, numScvIp,nodesPerElem,nDim);
    numScalars += numScvIp * nodesPerElem * nDim;
  }

  if (needDerivFem) {
    deriv_fem = get_shmem_view_3D<T>(team, numFemIp,nodesPerElem,nDim);
    numScalars += numFemIp * nodesPerElem * nDim;
  }

  if (needDetjFC) {
    det_j_fc_scs = get_shmem_view_1D<T>(team, numFaceIp);
    numScalars += numFaceIp;
  }

  if (needDetj) {
    det_j = get_shmem_view_1D<T>(team, numScsIp);
    numScalars += numScsIp;
  }

  if (needDetjScv) {
    det_j_scv = get_shmem_view_1D<T>(team, numScvIp);
    numScalars += numScvIp;
  }

  if (needDetjFem) {
    det_j_fem = get_shmem_view_1D<T>(team, numFemIp);
    numScalars += numFemIp;
  }

  // error check
  if ( femGradOp && femShiftedGradOp )
    ThrowRequireMsg(numFemIp > 0, "ERROR, femGradOp and femShiftedGradOp both requested.");

  return numScalars;
}

template<typename T>
void MasterElementViews<T>::fill_master_element_views(
  const std::set<ELEM_DATA_NEEDED>& dataEnums,
  SharedMemView<double**>* coordsView,
  MasterElement* meFC,
  MasterElement* meSCS,
  MasterElement* meSCV,
  MasterElement* meFEM,
  const int* faceOrdinals)
{
  // Guard against calling MasterElement methods on SIMD data structures
  static_assert(std::is_same<T, double>::value,
                "Cannot populate MasterElement Views for non-double data types");

  double error = 0.0;
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case FC_AREAV:
        ThrowRequireMsg(meFC != nullptr, "ERROR, meFC needs to be non-null if FC_AREAV is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but FC_AREAV requested.");
        meFC->determinant(1, &((*coordsView)(0,0)), &fc_areav(0,0), &error);
        break;
      case SCS_AREAV:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_AREAV is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_AREAV requested.");
        meSCS->determinant(1, &((*coordsView)(0, 0)), &scs_areav(0, 0), &error);
        break;
      case SCS_FACE_GRAD_OP:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_FACE_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        //faceOrdinals is an array of length stk::simd::ndoubles. For non-kokkos master-elem, just pick off first one:
        meSCS->face_grad_op(1, faceOrdinals[0], &((*coordsView)(0,0)), &dndx_fc_scs(0,0,0), &deriv_fc_scs(0,0,0), &det_j_fc_scs(0));
        break;
      case SCS_GRAD_OP:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        meSCS->grad_op(1, &((*coordsView)(0, 0)), &dndx(0, 0, 0), &deriv(0, 0, 0), &det_j(0), &error);
        break;
      case SCS_SHIFTED_GRAD_OP:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        meSCS->shifted_grad_op(1, &((*coordsView)(0, 0)), &dndx_shifted(0, 0, 0), &deriv(0, 0, 0), &det_j(0), &error);
        break;
      case SCS_GIJ:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GIJ is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GIJ requested.");
        meSCS->gij(&((*coordsView)(0, 0)), &gijUpper(0, 0, 0), &gijLower(0, 0, 0), &deriv(0, 0, 0));
        break;
      case SCV_VOLUME:
        ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_VOLUME is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCV_VOLUME requested.");
        meSCV->determinant(1, &((*coordsView)(0, 0)), &scv_volume(0), &error);
        break;
      case SCV_GRAD_OP:
        ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCV_GRAD_OP requested.");
        meSCV->grad_op(1, &((*coordsView)(0, 0)), &dndx_scv(0, 0, 0), &deriv_scv(0, 0, 0), &det_j_scv(0), &error);
        break;
      case FEM_GRAD_OP:
        ThrowRequireMsg(meFEM != nullptr, "ERROR, meFEM needs to be non-null if FEM_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but FEM_GRAD_OP requested.");
        meFEM->grad_op(1, &((*coordsView)(0, 0)), &dndx_fem(0, 0, 0), &deriv_fem(0, 0, 0), &det_j_fem(0), &error);
        break;
      case FEM_SHIFTED_GRAD_OP:
        ThrowRequireMsg(meFEM != nullptr, "ERROR, meFEM needs to be non-null if FEM_SHIFTED_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but FEM_GRAD_OP requested.");
        meFEM->shifted_grad_op(1, &((*coordsView)(0, 0)), &dndx_fem(0, 0, 0), &deriv_fem(0, 0, 0), &det_j_fem(0), &error);
        break;

      default:
        break;
    }
  }
}

template<typename T>
void MasterElementViews<T>::fill_master_element_views_new_me(
  const std::set<ELEM_DATA_NEEDED>& dataEnums,
  SharedMemView<DoubleType**>* coordsView,
  MasterElement* meFC,
  MasterElement* meSCS,
  MasterElement* meSCV,
  MasterElement* meFEM,
  const int* faceOrdinals)
{
  for(ELEM_DATA_NEEDED data : dataEnums) {
    switch(data)
    {
      case FC_AREAV:
        ThrowRequireMsg(false, "FC_AREAV not implemented yet.");
        break;
      case SCS_AREAV:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_AREAV is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_AREAV requested.");
         meSCS->determinant(*coordsView, scs_areav);
         break;
      case SCS_FACE_GRAD_OP:
        ThrowRequireMsg(false, "SCS_FACE_GRAD_OP not implemented yet.");
       break;
      case SCS_GRAD_OP:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
         meSCS->grad_op(*coordsView, dndx, deriv);
         break;
      case SCS_SHIFTED_GRAD_OP:
        ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        meSCS->shifted_grad_op(*coordsView, dndx_shifted, deriv);
        break;
      case SCS_GIJ:
         ThrowRequireMsg(meSCS != nullptr, "ERROR, meSCS needs to be non-null if SCS_GIJ is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GIJ requested.");
         meSCS->gij(*coordsView, gijUpper, gijLower, deriv);
         break;
      case SCV_VOLUME:
         ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_VOLUME is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCV_VOLUME requested.");
         meSCV->determinant(*coordsView, scv_volume);
         break;
      case SCV_GRAD_OP:
        ThrowRequireMsg(meSCV != nullptr, "ERROR, meSCV needs to be non-null if SCV_GRAD_OP is requested.");
        ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but SCS_GRAD_OP requested.");
        meSCV->grad_op(*coordsView, dndx_scv, deriv_scv);
        break;
      case FEM_GRAD_OP:
         ThrowRequireMsg(meFEM != nullptr, "ERROR, meFEM needs to be non-null if FEM_GRAD_OP is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but FEM_GRAD_OP requested.");
         meFEM->grad_op_fem(*coordsView, dndx_fem, deriv_fem, det_j_fem);
         break;
      case FEM_SHIFTED_GRAD_OP:
         ThrowRequireMsg(meFEM != nullptr, "ERROR, meFEM needs to be non-null if FEM_SHIFTED_GRAD_OP is requested.");
         ThrowRequireMsg(coordsView != nullptr, "ERROR, coords null but FEM_GRAD_OP requested.");
         meFEM->shifted_grad_op_fem(*coordsView, dndx_fem, deriv_fem, det_j_fem);
         break;

      default: break;
    }
  }
}

template<typename T>
ScratchViews<T>::ScratchViews(const TeamHandleType& team,
             const stk::mesh::BulkData& bulkData,
             int nodesPerEntity,
             const ElemDataRequests& dataNeeded)
{
  /* master elements are allowed to be null if they are not required */
  MasterElement *meFC = dataNeeded.get_cvfem_face_me();
  MasterElement *meSCS = dataNeeded.get_cvfem_surface_me();
  MasterElement *meSCV = dataNeeded.get_cvfem_volume_me();
  MasterElement *meFEM = dataNeeded.get_fem_volume_me();

  int nDim = bulkData.mesh_meta_data().spatial_dimension();
  int nodesPerFace = meFC != nullptr ? meFC->nodesPerElement_ : 0;
  int nodesPerElem = meSCS != nullptr
          ? meSCS->nodesPerElement_ : meSCV != nullptr
          ? meSCV->nodesPerElement_ : meFEM != nullptr
          ? meFEM->nodesPerElement_ : 0;
  int numFaceIp= meFC  != nullptr ? meFC->numIntPoints_  : 0;
  int numScsIp = meSCS != nullptr ? meSCS->numIntPoints_ : 0;
  int numScvIp = meSCV != nullptr ? meSCV->numIntPoints_ : 0;
  int numFemIp = meFEM != nullptr ? meFEM->numIntPoints_ : 0;

  create_needed_field_views(team, dataNeeded, bulkData, nodesPerEntity);

  create_needed_master_element_views(team, dataNeeded, nDim, nodesPerFace, nodesPerElem, numFaceIp, numScsIp, numScvIp, numFemIp);
}

template<typename T>
void ScratchViews<T>::create_needed_field_views(const TeamHandleType& team,
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
        fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<T*>>(get_shmem_view_1D<T>(team, scalarsDim1), 1);
        numScalars += scalarsDim1;
      }
      else {
        fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<T**>>(get_shmem_view_2D<T>(team, scalarsDim1, scalarsDim2),2);
        numScalars += scalarsDim1 * scalarsDim2;
      }
    }
    else if (fieldEntityRank==stk::topology::NODE_RANK) {
      if (scalarsDim2 == 0) {
        if (scalarsDim1 == 1) {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<T*>>(get_shmem_view_1D<T>(team, nodesPerElem),1);
          numScalars += nodesPerElem;
        }
        else {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<T**>>(get_shmem_view_2D<T>(team, nodesPerElem, scalarsDim1),2);
          numScalars += nodesPerElem * scalarsDim1;
        }
      }
      else {
          fieldViews[fieldInfo.field->mesh_meta_data_ordinal()] = new ViewT<SharedMemView<T***>>(get_shmem_view_3D<T>(team, nodesPerElem, scalarsDim1, scalarsDim2),3);
          numScalars += nodesPerElem * scalarsDim1 * scalarsDim2;
      }
    }
    else {
      ThrowRequireMsg(false,"Only elem-rank and node-rank fields supported for scratch-views currently.");
    }
  }

  // Track total bytes required for field allocations
  num_bytes_required += numScalars * sizeof(T);
}

template<typename T>
void ScratchViews<T>::create_needed_master_element_views(const TeamHandleType& team,
                                        const ElemDataRequests& dataNeeded,
                                        int nDim, int nodesPerFace, int nodesPerElem,
                                        int numFaceIp, int numScsIp, int numScvIp, int numFemIp)
{
  int numScalars = 0;

  for (auto it = dataNeeded.get_coordinates_map().begin();
       it != dataNeeded.get_coordinates_map().end(); ++it) {
    hasCoordField[it->first] = true;
    numScalars += meViews[it->first].create_master_element_views(
      team, dataNeeded.get_data_enums(it->first),
      nDim, nodesPerFace, nodesPerElem, numFaceIp, numScsIp, numScvIp, numFemIp);
  }

  num_bytes_required += numScalars * sizeof(T);
}

int get_num_scalars_pre_req_data( ElemDataRequests& dataNeededBySuppAlgs, int nDim);

void fill_pre_req_data(ElemDataRequests& dataNeeded,
                       const stk::mesh::BulkData& bulkData,
                       stk::mesh::Entity elem,
                       ScratchViews<double>& prereqData,
                       bool fillMEViews = true);

void fill_master_element_views(ElemDataRequests& dataNeeded,
                               const stk::mesh::BulkData& bulkData,
                               ScratchViews<DoubleType>& prereqData,
                               const int* faceOrdinals = nullptr);

template<typename T = double>
int get_num_bytes_pre_req_data(ElemDataRequests& dataNeededBySuppAlgs, int nDim)
{
  return sizeof(T) * get_num_scalars_pre_req_data(dataNeededBySuppAlgs, nDim);
}

} // namespace nalu
} // namespace Sierra

#endif
