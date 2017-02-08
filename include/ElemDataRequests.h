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
  NODES = 0,
  SCS_AREAV,
  SCS_GRAD_OP,
  SCV_VOLUME
};

struct FieldBaseLess
{
  bool operator()(const stk::mesh::FieldBase* lhs, const stk::mesh::FieldBase* rhs) const
  {
    return lhs->mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal();
  }
};

typedef std::set<const stk::mesh::FieldBase*,FieldBaseLess> FieldSet;

class ElemDataRequests
{
public:
  ElemDataRequests()
    : dataEnums(), fields(), meSCS_(NULL), meSCV_(NULL)
  {
  }

  void add(ELEM_DATA_NEEDED data)
  {
    dataEnums.insert(data);
  }
  void add_master_element_call(ELEM_DATA_NEEDED data)
  {
    dataEnums.insert(data);
  }

  void add(const stk::mesh::FieldBase& field)
  {
    fields.insert(&field);
  }
  void add_gathered_nodal_field(const stk::mesh::FieldBase& field)
  {
    fields.insert(&field);
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
