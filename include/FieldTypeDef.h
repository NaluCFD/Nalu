/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef FieldTypeDef_h
#define FieldTypeDef_h

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace sierra{
namespace nalu{

// define scalar field typedef
typedef stk::mesh::Field<double>  ScalarFieldType;
typedef stk::mesh::Field<stk::mesh::EntityId> GlobalIdFieldType;
typedef stk::mesh::Field<int>  ScalarIntFieldType;

// define vector field typedef; however, what is the value of Cartesian?
typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorFieldType;

// define generic
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;

// field type for local ids
typedef unsigned LocalId;
typedef stk::mesh::Field<LocalId>  LocalIdFieldType;

} // namespace nalu
} // namespace Sierra

#endif
