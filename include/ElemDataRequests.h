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
  SCS_AREAV = 0,
  SCS_GRAD_OP,
  SCV_VOLUME
};

struct FieldInfo {
  FieldInfo(const stk::mesh::FieldBase* fld, unsigned scalars)
  : field(fld), scalarsPerEntity(scalars)
  {}
  const stk::mesh::FieldBase* field;
  unsigned scalarsPerEntity;
};

struct FieldInfoLess
{
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
    : dataEnums(), fields(), meSCS_(NULL), meSCV_(NULL)
  {
  }

  void add_master_element_call(ELEM_DATA_NEEDED data)
  {
    dataEnums.insert(data);
  }

  void add_gathered_nodal_field(const stk::mesh::FieldBase& field, unsigned scalarsPerNode)
  {
    FieldInfo fieldInfo(&field, scalarsPerNode);
    FieldSet::iterator iter = fields.find(fieldInfo);
    if (iter == fields.end()) {
      fields.insert(FieldInfo(&field, scalarsPerNode));
    }
    else {
      ThrowAssertMsg(iter->scalarsPerEntity == scalarsPerNode, "ElemDataRequests ERROR, gathered-nodal-field "<<field.name()<<" requested with scalarsPerNode=="<<scalarsPerNode<<", but previously requested with scalarsPerNode=="<<iter->scalarsPerEntity);
    }
  }

  void add_cvfem_volume_me(MasterElement *meSCV)
  {
    meSCV_ = meSCV;
  }
  void add_cvfem_surface_me(MasterElement *meSCS)
  {
    meSCS_ = meSCS;
  }

  const std::set<ELEM_DATA_NEEDED>& get_data_enums() const { return dataEnums; }

  const FieldSet& get_fields() const { return fields; }  
  MasterElement *get_cvfem_volume_me() {return meSCV_;}
  MasterElement *get_cvfem_surface_me() {return meSCS_;}

private:
  std::set<ELEM_DATA_NEEDED> dataEnums;
  FieldSet fields;
  MasterElement *meSCS_;
  MasterElement *meSCV_;
};

} // namespace nalu
} // namespace Sierra

#endif
