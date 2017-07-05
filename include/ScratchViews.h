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

namespace sierra{
namespace nalu{

struct ViewHolder {
  virtual ~ViewHolder() {}
};

template<typename T>
struct ViewT : public ViewHolder {
  ViewT(T view) : view_(view) {}
  virtual ~ViewT(){}
  T view_;
};

class MasterElementViews
{
public:
  MasterElementViews() = default;
  virtual ~MasterElementViews() = default;

  int create_master_element_views(
    const TeamHandleType& team,
    const std::set<ELEM_DATA_NEEDED>& dataEnums,
    int nDim, int nodesPerElem,
    int numScsIp, int numScvIp, int numFemIp);

  void fill_master_element_views(
    const std::set<ELEM_DATA_NEEDED>& dataEnums,
    SharedMemView<double**>* coordsView,
    MasterElement* meSCS,
    MasterElement* meSCV,
    MasterElement* meFEM);

  SharedMemView<double**> scs_areav;
  SharedMemView<double***> dndx;
  SharedMemView<double***> dndx_shifted;
  SharedMemView<double***> dndx_fem;
  SharedMemView<double*> deriv;
  SharedMemView<double*> deriv_fem;
  SharedMemView<double*> det_j;
  SharedMemView<double*> det_j_fem;
  SharedMemView<double*> scv_volume;
  SharedMemView<double***> gijUpper;
  SharedMemView<double***> gijLower;
};

class ScratchViews
{
public:
  ScratchViews(const TeamHandleType& team,
               const stk::mesh::BulkData& bulkData,
               stk::topology topo,
               ElemDataRequests& dataNeeded);

  virtual ~ScratchViews() {
    for(ViewHolder* vh : fieldViews) {
      delete vh;
    }
  }

  inline
  SharedMemView<double*>& get_scratch_view_1D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<double**>& get_scratch_view_2D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<double***>& get_scratch_view_3D(const stk::mesh::FieldBase& field);

  inline
  SharedMemView<double****>& get_scratch_view_4D(const stk::mesh::FieldBase& field);

  inline
  MasterElementViews& get_me_views(const COORDS_TYPES cType)
  {
    ThrowRequire(hasCoordField[cType] == true);
    return meViews[cType];
  }

  inline int total_bytes() const { return num_bytes_required; }

  const stk::mesh::Entity* elemNodes;

private:
  void create_needed_field_views(const TeamHandleType& team,
                                 const ElemDataRequests& dataNeeded,
                                 const stk::mesh::BulkData& bulkData,
                                 int nodesPerElem);

  void create_needed_master_element_views(const TeamHandleType& team,
                                          const ElemDataRequests& dataNeeded,
                                          int nDim, int nodesPerElem,
                                          int numScsIp, int numScvIp, int numFemIp);

  std::vector<ViewHolder*> fieldViews;
  MasterElementViews meViews[MAX_COORDS_TYPES];
  bool hasCoordField[MAX_COORDS_TYPES] = {false, false};
  int num_bytes_required{0};
};

inline
SharedMemView<double*>& ScratchViews::get_scratch_view_1D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 1D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<double*>>* vt = static_cast<ViewT<SharedMemView<double*>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

inline
SharedMemView<double**>& ScratchViews::get_scratch_view_2D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 2D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<double**>>* vt = static_cast<ViewT<SharedMemView<double**>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

inline
SharedMemView<double***>& ScratchViews::get_scratch_view_3D(const stk::mesh::FieldBase& field)
{ 
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 3D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<double***>>* vt = static_cast<ViewT<SharedMemView<double***>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

inline
SharedMemView<double****>& ScratchViews::get_scratch_view_4D(const stk::mesh::FieldBase& field)
{
  ThrowAssertMsg(fieldViews[field.mesh_meta_data_ordinal()] != nullptr, "ScratchViews ERROR, trying to get 4D scratch-view for field "<<field.name()<<" which wasn't declared as pre-req field.");
  ViewT<SharedMemView<double****>>* vt = static_cast<ViewT<SharedMemView<double****>>*>(fieldViews[field.mesh_meta_data_ordinal()]);
  return vt->view_;
}

int get_num_bytes_pre_req_data( ElemDataRequests& dataNeededBySuppAlgs, int nDim);

void fill_pre_req_data(ElemDataRequests& dataNeeded,
                       const stk::mesh::BulkData& bulkData,
                       stk::topology topo,
                       stk::mesh::Entity elem,
                       ScratchViews& prereqData);

} // namespace nalu
} // namespace Sierra

#endif
