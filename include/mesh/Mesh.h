/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Mesh_h
#define Mesh_h

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>

#include<vector>

namespace stk {
namespace mesh {
class Part;
class Selector;
}
}

/** Nalu interface to STK mesh */

namespace sierra {
namespace nalu {
namespace mesh {

// put a field, of size 1, on a given part without a provided initial value
template< class field_type >
field_type &put_field( 
  field_type &field,
  const stk::mesh::Part &part) {
  const std::vector<typename stk::mesh::FieldTraits<field_type>::data_type> init_value(1,0);
  return stk::mesh::put_field_on_mesh(field, part, init_value.data());
 }

// put a field, of size n1, on a given part without a provided initial value
template< class field_type >
field_type & put_field( 
  field_type &field,
  const stk::mesh::Part &part,
  unsigned n1) {
  const std::vector<typename stk::mesh::FieldTraits<field_type>::data_type> init_value(n1,0);
  return stk::mesh::put_field_on_mesh(field, part, n1, init_value.data());
 }

// put a field, of size n1*n2, on a given part without a provided initial value
template< class field_type >
field_type & put_field( 
  field_type &field,
  const stk::mesh::Part &part,
  unsigned n1,
  unsigned n2) {
  const std::vector<typename stk::mesh::FieldTraits<field_type>::data_type> init_value(n1*n2,0);
  return stk::mesh::put_field_on_mesh(field, part, n1, n2, init_value.data());
 }

// put a field, of size 1, on a given part with a provided initial value
template< class field_type >
field_type &put_field_with_ic_value( 
  field_type &field,
  const stk::mesh::Part &part,
  const typename stk::mesh::FieldTraits<field_type>::data_type* init_value) {
  return stk::mesh::put_field_on_mesh(field, part, init_value);
 }

// put a field, of size n1, on a given part with a provided initial value
template< class field_type >
field_type & put_field_with_ic_value( 
  field_type &field,
  const stk::mesh::Part &part,
  unsigned n1,
  const typename stk::mesh::FieldTraits<field_type>::data_type* init_value) {
  return stk::mesh::put_field_on_mesh(field, part, n1, init_value);
 }

// put a field, of size n1*n2, on a given part with a provided initial value
template< class field_type >
field_type & put_field_with_ic_value( 
  field_type &field,
  const stk::mesh::Part &part,
  unsigned n1,
  unsigned n2,
  const typename stk::mesh::FieldTraits<field_type>::data_type* init_value) {
  return stk::mesh::put_field_on_mesh(field, part, n1, n2, init_value);
 }

// put a field, of size 1, on a given selector without a provided initial value
template< class field_type >
field_type &put_field( 
  field_type &field,
  const stk::mesh::Selector &selector) {
  const std::vector<typename stk::mesh::FieldTraits<field_type>::data_type> init_value(1,0);
  return stk::mesh::put_field_on_mesh(field, selector, init_value.data());
 }

// put a field, of size n1, on a given selector without a provided initial value
template< class field_type >
field_type & put_field( 
  field_type &field,
  const stk::mesh::Selector &selector,
  unsigned n1) {
  const std::vector<typename stk::mesh::FieldTraits<field_type>::data_type> init_value(n1,0);
  return stk::mesh::put_field_on_mesh(field, selector, n1, init_value.data());
 }

// put a field, of size 1, on a given selector with a provided initial value
template< class field_type >
field_type &put_field_with_ic_value( 
  field_type &field,
  const stk::mesh::Selector &selector,
  const typename stk::mesh::FieldTraits<field_type>::data_type* init_value) {
  return stk::mesh::put_field_on_mesh(field, selector, init_value);
 }

// put a field, of size n1, on a given selector with a provided initial value
template< class field_type >
field_type & put_field_with_ic_value( 
  field_type &field,
  const stk::mesh::Selector &selector,
  unsigned n1,
  const typename stk::mesh::FieldTraits<field_type>::data_type* init_value) {
  return stk::mesh::put_field_on_mesh(field, selector, n1, init_value);
 }

} // namespace mesh
} // namespace nalu
} // namespace Sierra

#endif
