/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ElemDataRequests_h
#define ElemDataRequests_h

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <set>

namespace sierra{
namespace nalu{

class MasterElement;

enum ELEM_DATA_NEEDED {
  FC_AREAV = 0,
  SCS_AREAV,
  SCS_FACE_GRAD_OP,
  SCS_SHIFTED_FACE_GRAD_OP,
  SCS_GRAD_OP,
  SCS_SHIFTED_GRAD_OP,
  SCS_GIJ,
  SCV_VOLUME,
  SCV_GRAD_OP,
  SCV_SHIFTED_GRAD_OP,
  FEM_GRAD_OP,
  FEM_SHIFTED_GRAD_OP,
  FEM_DET_J,
  FEM_NORMAL,
  FEM_FACE_NORMAL,
  FEM_FACE_GRAD_OP,
  FEM_FACE_DET_J,
  FEM_GIJ,
};

enum COORDS_TYPES {
  CURRENT_COORDINATES = 0,
  MODEL_COORDINATES,
  MAX_COORDS_TYPES,
};

static const std::string CoordinatesTypeNames[] = {
  "current_coordinates",
  "model_coordinates"
};

struct FieldInfo {
  FieldInfo(const stk::mesh::FieldBase* fld, unsigned scalars)
  : field(fld), scalarsDim1(scalars), scalarsDim2(0)
  {}
  FieldInfo(const stk::mesh::FieldBase* fld, unsigned tensorDim1, unsigned tensorDim2)
  : field(fld), scalarsDim1(tensorDim1), scalarsDim2(tensorDim2)
  {}
  const stk::mesh::FieldBase* field;
  unsigned scalarsDim1;
  unsigned scalarsDim2;
};

struct FieldInfoLess {
  bool operator()(const FieldInfo& lhs, const FieldInfo& rhs) const
  {
    return lhs.field->mesh_meta_data_ordinal() < rhs.field->mesh_meta_data_ordinal();
  }
};

typedef std::set<FieldInfo,FieldInfoLess> FieldSet;

class ElemDataRequests
{
public:
  ElemDataRequests()
    : dataEnums(),
      coordsFields_(),
    fields(), meFC_(nullptr), meSCS_(nullptr), meSCV_(nullptr), meFCFEM_(nullptr), meFEM_(nullptr)
  {
  }

  void add_master_element_call(
    ELEM_DATA_NEEDED data,
    COORDS_TYPES cType = CURRENT_COORDINATES
  )
  {
    auto it = coordsFields_.find(cType);
    ThrowRequireMsg(
      it != coordsFields_.end(),
      "ElemDataRequests::add_master_element_call: Coordinates field (" +
      CoordinatesTypeNames[cType] +
      ") must be registered to ElemDataRequests before registering MasterElement calls");
    dataEnums[cType].insert(data);
  }

  void add_gathered_nodal_field(const stk::mesh::FieldBase& field, unsigned scalarsPerNode);

  void add_gathered_nodal_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2);

  void add_face_field(const stk::mesh::FieldBase& field, unsigned scalarsPerFace);
  void add_ip_field(const stk::mesh::FieldBase& field, unsigned scalarsPerElement);
  void add_element_field(const stk::mesh::FieldBase& field, unsigned scalarsPerElement);

  void add_face_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2);
  void add_ip_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2);
  void add_element_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2);

  void add_coordinates_field(
    const stk::mesh::FieldBase& field,
    unsigned scalarsPerNode,
    COORDS_TYPES cType);

  void add_cvfem_face_me(MasterElement *meFC)
  {
    meFC_ = meFC;
  }

  void add_cvfem_volume_me(MasterElement *meSCV)
  {
    meSCV_ = meSCV;
  }

  void add_cvfem_surface_me(MasterElement *meSCS)
  {
    meSCS_ = meSCS;
  }

  void add_fem_face_me(MasterElement *meFCFEM)
  {
    meFCFEM_ = meFCFEM;
  }

  void add_fem_volume_me(MasterElement *meFEM)
  {
    meFEM_ = meFEM;
  }

  const std::set<ELEM_DATA_NEEDED>& get_data_enums(
    const COORDS_TYPES cType) const
  { return dataEnums[cType]; }

  const stk::mesh::FieldBase* get_coordinates_field(
    const COORDS_TYPES cType) const
  {
    auto it = coordsFields_.find(cType);
    ThrowRequireMsg(
      it != coordsFields_.end(),
      "ElemDataRequests::get_coordinates_field: Coordinates field (" +
      CoordinatesTypeNames[cType] +
      ") must be registered to ElemDataRequests before access");
    return it->second;
  }

  const std::map<COORDS_TYPES, const stk::mesh::FieldBase*>&
  get_coordinates_map() const
  { return coordsFields_; }

  const FieldSet& get_fields() const { return fields; }  
  MasterElement *get_cvfem_face_me() const {return meFC_;}
  MasterElement *get_cvfem_volume_me() const {return meSCV_;}
  MasterElement *get_cvfem_surface_me() const {return meSCS_;}
  MasterElement *get_fem_volume_me() const {return meFEM_;}
  MasterElement *get_fem_face_me() const {return meFCFEM_;}

private:
  std::array<std::set<ELEM_DATA_NEEDED>, MAX_COORDS_TYPES> dataEnums;
  std::map<COORDS_TYPES, const stk::mesh::FieldBase*> coordsFields_;
  FieldSet fields;
  MasterElement *meFC_;
  MasterElement *meSCS_;
  MasterElement *meSCV_;
  MasterElement *meFCFEM_;
  MasterElement *meFEM_;
};

} // namespace nalu
} // namespace Sierra

#endif
