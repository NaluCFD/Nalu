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
  : field(fld), scalarsDim1(scalars), scalarsDim2(0)
  {}
  FieldInfo(const stk::mesh::FieldBase* fld, unsigned tensorDim1, unsigned tensorDim2)
  : field(fld), scalarsDim1(tensorDim1), scalarsDim2(tensorDim2)
  {}
  const stk::mesh::FieldBase* field;
  unsigned scalarsDim1;
  unsigned scalarsDim2;
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
    ThrowAssertMsg(field.entity_rank()==stk::topology::NODE_RANK,"ElemDataRequests ERROR, add_gathered_nodal_field called with field "<<field.name()<<" which has entity-rank=="<<field.entity_rank()<<" but is expected to have rank==NODE_RANK");

    FieldInfo fieldInfo(&field, scalarsPerNode);
    FieldSet::iterator iter = fields.find(fieldInfo);
    if (iter == fields.end()) {
      fields.insert(fieldInfo);
    }
    else {
      ThrowAssertMsg(iter->scalarsDim1 == scalarsPerNode, "ElemDataRequests ERROR, gathered-nodal-field "<<field.name()<<" requested with scalarsPerNode=="<<scalarsPerNode<<", but previously requested with scalarsPerNode=="<<iter->scalarsDim1);
    }
  }

  void add_gathered_nodal_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2)
  {
    ThrowAssertMsg(field.entity_rank()==stk::topology::NODE_RANK,"ElemDataRequests ERROR, add_gathered_nodal_field called with field "<<field.name()<<" which has entity-rank=="<<field.entity_rank()<<" but is expected to have rank==NODE_RANK");

    FieldInfo fieldInfo(&field, tensorDim1, tensorDim2);
    FieldSet::iterator iter = fields.find(fieldInfo);
    if (iter == fields.end()) {
      fields.insert(fieldInfo);
    }
    else {
      ThrowAssertMsg(iter->scalarsDim1 == tensorDim1 && iter->scalarsDim2 == tensorDim2, "ElemDataRequests ERROR, gathered-nodal-field "<<field.name()<<" requested with tensorDim1=="<<tensorDim1<<",tensorDim2=="<<tensorDim2<<", but previously requested with tensorDim1=="<<iter->scalarsDim1<<",tensorDim2=="<<iter->scalarsDim2);
    }
  }

  void add_element_field(const stk::mesh::FieldBase& field, unsigned scalarsPerElement)
  {
    ThrowAssertMsg(field.entity_rank()==stk::topology::ELEM_RANK,"ElemDataRequests ERROR, add_element_field called with field "<<field.name()<<" which has entity-rank=="<<field.entity_rank()<<" but is expected to have rank==ELEM_RANK");

    FieldInfo fieldInfo(&field, scalarsPerElement);
    FieldSet::iterator iter = fields.find(fieldInfo);
    if (iter == fields.end()) {
      fields.insert(fieldInfo);
    }
    else {
      ThrowAssertMsg(iter->scalarsDim1 == scalarsPerElement, "ElemDataRequests ERROR, element-field "<<field.name()<<" requested with scalarsPerElement=="<<scalarsPerElement<<", but previously requested with scalarsPerElement=="<<iter->scalarsDim1);
    }
  }

  void add_element_field(const stk::mesh::FieldBase& field, unsigned tensorDim1, unsigned tensorDim2)
  {
    ThrowAssertMsg(field.entity_rank()==stk::topology::ELEM_RANK,"ElemDataRequests ERROR, add_element_field called with field "<<field.name()<<" which has entity-rank=="<<field.entity_rank()<<" but is expected to have rank==ELEM_RANK");

    FieldInfo fieldInfo(&field, tensorDim1, tensorDim2);
    FieldSet::iterator iter = fields.find(fieldInfo);
    if (iter == fields.end()) {
      fields.insert(fieldInfo);
    }
    else {
      ThrowAssertMsg(iter->scalarsDim1 == tensorDim1 && iter->scalarsDim2 == tensorDim2, "ElemDataRequests ERROR, element-field "<<field.name()<<" requested with tensorDim1=="<<tensorDim1<<",tensorDim2=="<<tensorDim2<<", but previously requested with tensorDim1=="<<iter->scalarsDim1<<",tensorDim2=="<<iter->scalarsDim2);
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
